#!/usr/bin/env -S awk -f

# This awk script splits a .pheno file containing multiple phenotypes into a set of
# files, each with one of the phenotypes in it. The output files are written to the
# current directory

# param1: The path to a .pheno file containing the phenotype data

BEGIN {
    FS = OFS = "\t"
}

# Process the header to get phenotype names and prepare files
NR == 1 {
    for (i = 2; i <= NF; i++) {
        header[i] = $i;  # Store phenotype names
        file[i] = $i ".pheno";  # Construct filename for each phenotype
        print "#IID", $i > file[i]  # Print header to each file
    }
}

# Process the data rows
NR > 1 {
    for (i = 2; i <= NF; i++) {
        print $1, $i > file[i]  # Print IID and phenotype value to respective file
    }
}

END {
    for (i in file) {
        close(file[i]);  # Close all files
    }
}
