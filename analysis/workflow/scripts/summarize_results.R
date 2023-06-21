#!/usr/bin/env Rscript

# This R script summarizes the output of several executions of FINEMAP and SuSiE
# It cannot be executed independently of Snakemake






# handle error messages properly:
# they should all get written to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log, type='message')


suppressMessages(library(ggplot2))
suppressMessages(library(pgenlibr))
suppressMessages(library(data.table))

readPGEN = function(pfile) {
  if (endsWith(pfile, ".pgen")) {
    pfile = substr(pfile, 1, nchar(pfile)-5)
  }
  pvar = pgenlibr::NewPvar(paste0(pfile, '.pvar'))
  pgen = pgenlibr::NewPgen(paste0(pfile, '.pgen'), pvar=pvar)
  n = pgenlibr::GetVariantCt(pgen)
  X = pgenlibr::ReadList(pgen, 1:n, meanimpute=FALSE)
  pgenlibr::ClosePgen(pgen)
  pgenlibr::ClosePvar(pvar)
  X
}

write("Loading input data", stderr())
# import the finemap and susie results
# and the path to an output directory
gt = readPGEN(snakemake@input[["gt"]])
if ("finemap" %in% names(snakemake@input)) {
    finemap_results = readRDS(snakemake@input[["finemap"]])[[1]]
} else {
    finemap_results = NULL
}
susie_results = readRDS(snakemake@input[["susie"]])
out = snakemake@params[["outdir"]]
hap = FALSE
if (nchar(snakemake@params[["causal_hap"]]) > 0) {
    hap = TRUE
    happler_haplotype = file(snakemake@input[["happler_hap"]], "rt")
    causal_haplotype = file(snakemake@params[["causal_hap"]], "rt")
}

if (hap) {
    write("Parsing hap files", stderr())
    while(TRUE) {
        line = readLines(happler_haplotype, 1)
        if (grepl("^H", line)) break
    }
    # extract the start and end of the haplotype
    happler_hap = as.integer(strsplit(line, "\t")[[1]][c(3,4)])
    names(happler_hap) = c("start", "end")
    while(TRUE) {
        line = readLines(causal_haplotype, 1)
        if (grepl("^H", line)) break
    }
    # extract the start, end, and ID of the haplotype
    causal_hap = as.integer(strsplit(line, "\t")[[1]][c(3,4)])
    names(causal_hap) = c("start", "end")
    haplotypes = t(data.frame(causal_hap, happler_hap))
}

write("Formatting genotypes and creating output dir", stderr())
storage.mode(X) = 'double'
dir.create(out, showWarnings = FALSE)

# parse the susie data:
# 1) the truly causal variant, as defined in the simulation
# 2) whether the causal variant was provided in the genotypes
# 3) the list provided by susie as output
# 4) the PIPs output by SuSiE
causal_variant = susie_results$causal_var
exclude_causal = susie_results$causal_excluded
fitted = susie_results$fitted
susie_pip = fitted$pip

write("Handling causal variable", stderr())
# first, we must create a vector with the causal status
# this vector indicates which variant is truly causal
b = rep(0, ncol(X))
names(b) = colnames(X)
if (exclude_causal) {
    b = b[!(names(b) %in% c(causal_variant))]
    if (hap) {
        haplotypes = t(data.frame(haplotypes[c(2),]))
    }
} else {
    b[causal_variant] = 1
}

if (!is.null(finemap_results)) {
    write("Parsing FINEMAP data", stderr())
    # parse the finemap data
    snp = finemap_results$snp
    finemap_pip = snp[order(as.numeric(snp$snp)),]$snp_prob
    names(finemap_pip) = names(susie_pip)
}

# define a function that generates the PIP plot data
pip_plot_data = function(pips, X, b, susie_cs=NULL) {
    # note that each of pips, X, and b must be ordered by variant POS
    # first, initialize the values we need
    b_colors = c(`0`='black', `1`='red')
    causal_var = names(b[b == 1])
    if (length(causal_var) == 0) {
        causal_var = X[,causal_variant]
        X = X[,!(colnames(X) %in% c(causal_variant))]
    } else {
        causal_var = X[,causal_var[1]]
    }
    data = data.frame(
        pip = pips,
        b = as.character(b),
        pos = as.integer(sub("\\.(0|1|2)", "", sub("X(\\S+)", "\\1", names(b)))),
        ld_causal = as.vector(cor(causal_var, X))^2
    )
    if (!is.null(susie_cs)) {
        data$cs = as.integer((names(b) %in% susie_cs))*2
    } else {
        data$cs = 0
    }
    return(data)
}

# define a function for PIP plotting
pip_plot = function(pips, X, b, susie_cs=NULL) {
    # create a ggplot of the PIPs
    # but first, get the data we need for the plot
    data = pip_plot_data(pips, X, b, susie_cs)
    # extract the causal variants to another data frame
    data_causal = data[data$b == '1',]
    # make the plot
    ggplot(data, aes(x=pos, y=pip)) +
    geom_point(aes(fill=ld_causal, stroke=cs, color=factor(cs)), size=7, shape=21) +
    scale_fill_gradient(name='LD with Causal Variant', low='#FBBA72', high='#691E06') +
    scale_color_manual(name='Credible Sets', values=c('transparent', '#7C9299'), guide="none") +
    geom_point(data=data_causal, aes(stroke=cs, color=factor(cs)), fill='red', size=7, shape=21) +
    xlab('Chromosomal Position') +
    ylab('Posterior Inclusion Probability (PIP)') + 
    ylim(0,1) +
    theme_grey(base_size=16)
}

# define a function for PIP plotting with haplotypes as lines
pip_plot_haps = function(pips, X, b, haplotypes, susie_cs=NULL) {
    # create a ggplot of the PIPs
    # but first, get the data we need for the plot
    data = pip_plot_data(pips, X, b, susie_cs)
    # extract the haplotypes to another data frame
    data_hap = cbind(data[grepl(".2", rownames(data), fixed=T),], haplotypes)
    data_hap$color = c("black", "red")[as.integer(data_hap$b)+1]
    # remove the haps from the data
    data = data[!(row.names(data) %in% rownames(data_hap)),]
    # make the plot
    ggplot(data, aes(x=pos, y=pip)) +
    geom_point(aes(fill=ld_causal, stroke=cs, color=factor(cs)), size=7, shape=21) +
    scale_fill_gradient(name='LD with Causal Variant', low='#FBBA72', high='#691E06') +
    scale_color_manual(name='Credible Sets', values=c('transparent', '#7C9299'), guide="none") +
    geom_segment(data=data_hap, aes(x = start, xend = end, y = pip, yend = pip, color=factor(b)), color=data_hap$color, size=3) +
    xlab('Chromosomal Position') +
    ylab('Posterior Inclusion Probability (PIP)') +
    ylim(0,1) +
    theme_grey(base_size=16)
}

# # plot the results of susie
# pdf(paste0(out,'/susie.pdf'), width =5, height = 5, pointsize=16)
# susieR::susie_plot(fitted, y='PIP', b=b, max_cs=0.4, main = paste('SuSiE, ', length(fitted$sets$cs), 'CS identified'))
# dev.off()

# bhat = coef(fitted)[-1]
# pdf(paste0(out,'/susie_eff.pdf'), width =5, height = 5, pointsize=16)
# susieR::susie_plot(bhat, y='bhat', b=b, main = 'SuSiE, effect size estimate') 
# dev.off()

# # plot the results of running FINEMAP
# pdf(paste0(out,'/finemap.pdf'), width =5, height = 5, pointsize=16)
# susieR::susie_plot(pip, y='PIP', b=b, main = 'Bayesian sparse regression')
# dev.off()


write("Creating PIP plot for SuSiE", stderr())
if (hap) {
    pip_plot_haps(susie_pip, X, b, haplotypes, susie_cs=names(susie_pip[fitted$sets$cs[['L1']]]))
} else {
    pip_plot(susie_pip, X, b, susie_cs=names(susie_pip[fitted$sets$cs[['L1']]]))
}
ggsave(paste0(out,'/susie.pdf'), width=10, height=5, device='pdf')
dev.off()

if (!is.null(finemap_results)) {
    write("Creating PIP plot for FINEMAP", stderr())
    if (hap) {
        pip_plot_haps(finemap_pip, X, b, haplotypes)
    } else {
        pip_plot(finemap_pip, X, b)
    }
    ggsave(paste0(out,'/finemap.pdf'), width=10, height=5, device='pdf')
    while (!is.null(dev.list())) dev.off()
}
