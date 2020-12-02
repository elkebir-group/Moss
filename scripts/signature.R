library('deconstructSigs')

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2 && length(args) != 3) {
    stop("Usage: signature.R input.csv output.csv [hg38]", call.=FALSE)
}

ifile <- args[1]
ofile <- args[2]
is_hg38 <- FALSE
if (length(args) == 3 && args[3] == "hg38") {
    is_hg38 <- TRUE
    print("Using hg38 genome.")
} else {
    print("Using hg19 genome.")
}

callset.df <- read.csv(file=ifile, header=TRUE, sep=",")
if (!is_hg38) {
    callset.df$chr <- paste("chr", callset.df$chr, sep='')
}
callset.df$sample <- "None"

if (is_hg38) {
    callset.input <- mut.to.sigs.input(mut.ref = callset.df,
                                    sample.id = "sample",
                                    chr = "chr",
                                    pos = "pos",
                                    ref = "ref",
                                    alt = "alt")
} else {
    callset.input <- mut.to.sigs.input(mut.ref = callset.df,
                                sample.id = "sample",
                                chr = "chr",
                                pos = "pos",
                                ref = "ref",
                                alt = "alt")
}

if (is_hg38) {
    callset.output = whichSignatures(tumor.ref = callset.input,
                                    signatures.ref = signatures.cosmic,
                                    sample.id = 'None',
                                    contexts.needed = TRUE,
                                    tri.counts.method = 'genome',
                                    signature.cutoff = 0)
} else {
    callset.output = whichSignatures(tumor.ref = callset.input,
                                    signatures.ref = signatures.cosmic,
                                    sample.id = 'None',
                                    contexts.needed = TRUE,
                                    tri.counts.method = 'exome2genome',
                                    signature.cutoff = 0)
}

write.csv(callset.output$weights, file=ofile)