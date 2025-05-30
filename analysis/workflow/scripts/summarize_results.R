#!/usr/bin/env Rscript

# This R script summarizes the output of several executions of FINEMAP and SuSiE
# It cannot be executed independently of Snakemake






# handle error messages properly:
# they should all get written to the log file
log <- file(snakemake@log[[1]], open="wt")
sink(log, type='message')

suppressMessages(library(ggplot2))

# load functions to help read various file types
source("workflow/scripts/utils.R")

write("Loading input data", stderr())
region = snakemake@params[["region"]]
phen = snakemake@input[["phen"]]
samples = readPSAM(snakemake@input[["gt"]], samples=readPheno(phen)[,1])
# import the finemap and susie results
# and the path to an output directory
X = readPGEN(snakemake@input[["gt"]], region=region, samples=samples[,2])
# also load the positions of each of the variants
pos = readPVAR(snakemake@input[["gt"]], region=region)
if (length(snakemake@input[["causal_gt"]])) {
    write("Loading causal genotype data too", stderr())
    causal_gt_file = snakemake@input[["causal_gt"]]
    causal_gt_samples = readPSAM(causal_gt_file, samples=readPheno(phen)[,1])
    causal_gt = readPGEN(causal_gt_file, region=region, samples=causal_gt_samples[,2])
    X = cbind(X, causal_gt)
    causal_variant = colnames(causal_gt)[1]
    # also load the positions of each of the variants
    pos = c(pos, readPVAR(causal_gt_file, region=region))
    write(paste("Loaded", length(pos), "positions"), stderr())
} else {
    causal_variant = NULL
    write("Causal haplotype is unknown", stderr())
}

if ("finemap" %in% names(snakemake@input)) {
    finemap_results = readRDS(snakemake@input[["finemap"]])[[1]]
} else {
    finemap_results = NULL
}
susie_results = readRDS(snakemake@input[["susie"]])
out = snakemake@params[["outdir"]]
hap = FALSE
causal_hap_is_provided = FALSE
if (length(snakemake@input[["happler_hap"]]) > 0) {
    hap = TRUE
    happler_haplotype = snakemake@input[["happler_hap"]]
}
if (length(snakemake@params[["causal_hap"]]) > 0) {
    causal_hap_is_provided = TRUE
    causal_haplotype = snakemake@params[["causal_hap"]]
}

if (hap) {
    write("Parsing happler hap file", stderr())
    happler_hap = readHap(happler_haplotype)
}
if (causal_hap_is_provided) {
    write("Parsing causal hap file", stderr())
    causal_hap = readHap(causal_haplotype)
    stopifnot(causal_variant %in% causal_hap["id"])
    haplotypes = cbind(causal_hap, happler_hap)
} else if (hap) {
    haplotypes = happler_hap
}

write("Formatting genotypes and creating output dir", stderr())
storage.mode(X) = 'double'
dir.create(out, showWarnings = FALSE)

# parse the susie data:
# 1) whether the causal variant was provided in the genotypes
# 2) the list provided by susie as output
# 3) the PIPs output by SuSiE
exclude_causal = (susie_results$causal_excluded == "NULL")
fitted = susie_results$fitted
susie_pip = fitted$pip

write("Handling causal variable", stderr())
# first, we must create a vector with the causal status
# this vector indicates which variant is truly causal
b = rep(0, ncol(X))
names(b) = colnames(X)
if (!is.null(causal_variant)) {
    if (exclude_causal) {
        write("Excluding causal variable", stderr())
        causal_var_idx = which.max(names(b) %in% c(causal_variant))
        b = b[-causal_var_idx]
        pos = pos[-causal_var_idx]
        if (hap && !causal_hap_is_provided) {
            haplotypes = t(data.frame(haplotypes[c(2),]))
        }
        stopifnot(nrow(haplotypes) == 1)
    } else {
        b[causal_variant] = 1
    }
} else {
    write("Causal variable already excluded", stderr())
}

if (!is.null(finemap_results)) {
    write("Parsing FINEMAP data", stderr())
    # parse the finemap data
    snp = finemap_results$snp
    finemap_pip = snp[order(as.numeric(snp$snp)),]$snp_prob
    names(finemap_pip) = names(susie_pip)
}

# define a function that generates the PIP plot data
pip_plot_data = function(pips, X, b, pos, susie_cs=NULL) {
    write("Creating dataframe for PIP plot", stderr())
    # note that each of pips, X, and b must be ordered by variant POS
    # first, initialize the values we need
    b_colors = c(`0`='black', `1`='red')
    causal_var = names(b[b == 1])
    if (!is.null(causal_variant)) {
        if (length(causal_var) == 0) {
            causal_var = X[,causal_variant]
            X = X[,-which.max(colnames(X) %in% c(causal_variant))]
        } else {
            # TODO: make sure this only outputs a single row?
            causal_var = X[,causal_var[1]]
        }
        ld_causal = as.vector(cor(causal_var, X))^2
        stopifnot(length(b) == length(pos))
        stopifnot(names(b) == names(pos))
        stopifnot(names(b) == names(pips))
        data = data.frame(
            pip = pips,
            b = as.character(b),
            pos = pos,
            ld_causal = ld_causal
        )
    } else {
        data = data.frame(
            pip = pips,
            b = as.character(b),
            pos = pos
        )
    }
    if (!is.null(susie_cs)) {
        data$cs = as.integer((names(b) %in% susie_cs))*2
    } else {
        data$cs = 0
    }
    write("Finished dataframe for PIP plot", stderr())
    return(data)
}

# define a function for PIP plotting
pip_plot = function(pips, X, b, pos, susie_cs=NULL) {
    # create a ggplot of the PIPs
    # but first, get the data we need for the plot
    data = pip_plot_data(pips, X, b, pos, susie_cs)
    # extract the causal variants to another data frame
    data_causal = data[data$b == '1',]
    # make the plot
    plt = ggplot(data, aes(x=pos, y=pip))
    if ("ld_causal" %in% colnames(data)) {
        write("Plotting with causal LD", stderr())
        plt = plt + geom_point(aes(fill=ld_causal, stroke=cs, color=factor(cs)), size=7, shape=21) +
        scale_fill_gradient(name='LD with Causal Variant', low='#FBBA72', high='#691E06')
    } else {
        write("Not including causal variable", stderr())
        plt = plt + geom_point(aes(fill=1, stroke=cs, color=factor(cs)), size=7, shape=21, show.legend=F)
    }
    plt + scale_color_manual(name='Credible Sets', values=c('transparent', '#7C9299'), guide="none") +
    geom_point(data=data_causal, aes(stroke=cs, color=factor(cs)), fill='red', size=7, shape=21) +
    xlab('Chromosomal Position') +
    ylab('Posterior Inclusion Probability (PIP)') + 
    ylim(0,1) +
    theme_grey(base_size=16)
}

# define a function for PIP plotting with haplotypes as lines
pip_plot_haps = function(pips, X, b, pos, haplotypes, susie_cs=NULL) {
    # create a ggplot of the PIPs
    # but first, get the data we need for the plot
    data = pip_plot_data(pips, X, b, pos, susie_cs)
    # extract start and end coords
    # note that we must convert to integer b/c they're characters
    start_end = t(data.frame(as.integer(haplotypes[1, c("start", "end")])))
    colnames(start_end) = c("start", "end")
    # extract the haplotypes to another data frame
    data_hap = cbind(data[haplotypes[1, "id"],], start_end)
    data_hap$color = c("black", "red")[as.integer(data_hap$b)+1]
    # remove the haps from the data
    data = data[!(row.names(data) %in% rownames(data_hap)),]
    # make the plot
    plt = ggplot(data, aes(x=pos, y=pip))
    if ("ld_causal" %in% colnames(data)) {
        write("Plotting with causal LD", stderr())
        plt = plt + geom_point(aes(fill=ld_causal, stroke=cs, color=factor(cs)), size=7, shape=21) +
        scale_fill_gradient(name='LD with Causal Variant', low='#FBBA72', high='#691E06')
    } else {
        write("Not including causal variable", stderr())
        plt = plt + geom_point(aes(fill=1, stroke=cs, color=factor(cs)), size=7, shape=21, show.legend=F)
    }
    plt + scale_color_manual(name='Credible Sets', values=c('transparent', '#7C9299'), guide="none") +
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
    pip_plot_haps(susie_pip, X, b, pos, haplotypes, susie_cs=names(susie_pip[fitted$sets$cs[['L1']]]))
} else {
    pip_plot(susie_pip, X, b, pos, susie_cs=names(susie_pip[fitted$sets$cs[['L1']]]))
}
ggsave(paste0(out,'/susie.pdf'), width=6, height=5, device='pdf')
dev.off()

if (!is.null(finemap_results)) {
    write("Creating PIP plot for FINEMAP", stderr())
    if (hap) {
        pip_plot_haps(finemap_pip, X, b, pos, haplotypes)
    } else {
        pip_plot(finemap_pip, X, b, pos)
    }
    ggsave(paste0(out,'/finemap.pdf'), width=6, height=5, device='pdf')
    while (!is.null(dev.list())) dev.off()
}
