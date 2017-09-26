suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  make_option("--pdf", type="character", help="path and filename of the output pdf file"),
  make_option("--file", type="character", help="path and name of input file")
)

opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS]",
                           option_list = option_list,
                           add_help_option=TRUE)
opt <- parse_args(opt_parser)

DEBUG <- FALSE
#DEBUG <- TRUE
if(DEBUG == TRUE){
  opt = list()
  opt$file = "TIN_table.tsv"
  opt$pdf = "test.pdf"
}

input <- opt$file
out_pdf <- opt$pdf
nr_of_samples <- opt$samples
ylabel <- opt$ylab
main_text <- opt$main

###############################################################################


data <- read.table(input, header =T, row.names = 1)
data[data == 0.0] <- NA

num_of_cols <- dim(data)[2]
width_size <- num_of_cols * 0.35

if( width_size < 2 ) {
width_size <- 2
}

pdf(out_pdf, height=8, width= width_size)
par(mar=c(15.1,4.1, 4.1, 1.0))
boxplot(data, las=3, cex.axis=0.6)
dev.off()
