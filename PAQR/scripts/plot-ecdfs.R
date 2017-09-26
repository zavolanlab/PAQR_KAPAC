suppressPackageStartupMessages(library("optparse"))

################################################################################
# define and load arguments/params
################################################################################

option_list <- list(
  make_option(c("-f", "--file"), type="character", help="name of the input file"),
  make_option("--main", type="character", default="ECDF(x)", help="title of the plot [default: %default]"),
  make_option("--xlab", type="character", default="average exon length", help="x-axis label of the plot [default: %default]"),
  make_option("--ylab", type="character", default="cumulative fraction of exons", help="x-axis label of the plot [default: %default]"),
  make_option("--pdf", type="character", help="path and filename of the output pdf")
)

opt_parser <- OptionParser(usage="Usage: %prog [OPTIONS]",
                           option_list = option_list,
                           add_help_option=TRUE)
opt <- parse_args(opt_parser)

#debugging <- TRUE
debugging <- FALSE
if (debugging == TRUE)
{
  opt = list()
  opt$file = "fusessh_schmiral/collaborations/tcell_differentiation/05_t_cell_paper_test/relLengthPerExon.rep1.1minTPM.all.tsv"
  opt$main="ECDF(x)"
  opt$xlab="average exon length"
  opt$ylab="cumulative fraction of exons"
  opt$pdf = "test.pdf"
}

################################################################################
# define a function to get distinct colors
################################################################################

iwanthue <- function(n, hmin=0, hmax=360, cmin=0, cmax=180, lmin=0, lmax=100, 
                     plot=FALSE, random=FALSE) {
  # Presently doesn't allow hmax > hmin (H is circular)
  # n: number of colours
  # hmin: lower bound of hue (0-360)
  # hmax: upper bound of hue (0-360)
  # cmin: lower bound of chroma (0-180)
  # cmax: upper bound of chroma (0-180)
  # lmin: lower bound of luminance (0-100)
  # lmax: upper bound of luminance (0-100)
  # plot: plot a colour swatch?
  # random: should clustering be random? (if FALSE, seed will be set to 1,
  #         and the RNG state will be restored on exit.) 
  require(colorspace)
  stopifnot(hmin >= 0, cmin >= 0, lmin >= 0, 
            hmax <= 360, cmax <= 180, lmax <= 100, 
            hmin <= hmax, cmin <= cmax, lmin <= lmax,
            n > 0)
  if(!random) {
    if (exists(".Random.seed", .GlobalEnv)) {
      old_seed <- .GlobalEnv$.Random.seed
      on.exit(.GlobalEnv$.Random.seed <- old_seed)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv))
    }
    set.seed(1)
  }
  lab <- LAB(as.matrix(expand.grid(seq(0, 100, 1), 
                                   seq(-100, 100, 5), 
                                   seq(-110, 100, 5))))
  if (any((hmin != 0 || cmin != 0 || lmin != 0 ||
           hmax != 360 || cmax != 180 || lmax != 100))) {
    hcl <- as(lab, 'polarLUV')
    hcl_coords <- coords(hcl)
    hcl <- hcl[which(hcl_coords[, 'H'] <= hmax & hcl_coords[, 'H'] >= hmin &
                       hcl_coords[, 'C'] <= cmax & hcl_coords[, 'C'] >= cmin & 
                       hcl_coords[, 'L'] <= lmax & hcl_coords[, 'L'] >= lmin), ]
    #hcl <- hcl[-which(is.na(coords(hcl)[, 2]))]
    lab <- as(hcl, 'LAB')    
  }
  lab <- lab[which(!is.na(hex(lab))), ]
  clus <- kmeans(coords(lab), n, iter.max=50)
  if (isTRUE(plot)) {
    swatch(hex(LAB(clus$centers)))
  }
  hex(LAB(clus$centers))
}

################################################################################
# define function to infer the patient and the sample type from TCGA barcodes
################################################################################

get_tissue_code <- function(barcode){
substring(tail(strsplit(strsplit(barcode,"_")[[1]][2], "-")[[1]], n=1),1,2)
}

get_patient_id <- function(barcode){
head(tail(strsplit(strsplit(barcode,"_")[[1]][2], "-")[[1]], n=2),n=1)
}


################################################################################
################################################################################

input <- opt$file
plotMain <- opt$main
xlabel <- opt$xlab
ylabel <- opt$ylab
pdf_out <- opt$pdf

# load input table
lengthPerExon <- read.table(input, head=T, row.names = 1, check.names = FALSE)
# sort them by patient ids (sample from the same patient are next to each other)
patient_ids <- sapply(colnames(lengthPerExon), get_patient_id)
lengthPerExon <- lengthPerExon[, order(patient_ids)]
patient_ids <- patient_ids[order(patient_ids)]

# colnames contains the TCGA-barcodes
# check all cols that are from primary tumor samples

tissue_code <- sapply(colnames(lengthPerExon), get_tissue_code)

################################################################################
# define the set of used colors and the color for each column
################################################################################

# the number of needed colors equals the number of distinct patients
nr_colors <- length(unique( patient_ids))

# if( nr_colors >= 2) {
# colors <- rainbow(nr_colors)
# } else {
# colors <- rainbow(2)
# }

# use johnbaums' code of iwanthue to get distinct colors
colors <- iwanthue(nr_colors)

# preset the color per column as unset
colors_per_columns <- rep("unset", dim(lengthPerExon)[2])
# lty_per_columns <- rep(2,dim(lengthPerExon)[2])

# define the color for each column (based on the patient id)
tmp_cnt <- 0
patients_visited_already <- c()
for (n in 1:dim(lengthPerExon)[2]) {
if( patient_ids[n] %in% patients_visited_already ) {
    colors_per_columns[n] <- colors[tmp_cnt]
    colors_per_columns[n] <- colors[ which(patients_visited_already == patient_ids[n]) ]
} else {
  tmp_cnt <- tmp_cnt + 1
  colors_per_columns[n] <- colors[tmp_cnt]
  patients_visited_already[tmp_cnt] <- patient_ids[n]
}

# if(colors_per_columns[n] != "unset") next
# tmp_cnt <- tmp_cnt + 1
# curr_col_name <- colnames(lengthPerExon)[n]
# curr_type <- tissue_code[n]
# patient <- strsplit(curr_col_name,"-")[[1]][3]
# tmp_pair <- grep( patient, colnames(lengthPerExon) )
# stopifnot(length(tmp_pair) == 2)
# matching_sample_idx <- tmp_pair[tmp_pair != n]
# colors_per_columns[n] <- colors[tmp_cnt]
# colors_per_columns[matching_sample_idx] <- colors[tmp_cnt]
# if(grepl("^01",curr_type)){
# lty_per_columns[ matching_sample_idx ] <- 2
# } else{
# lty_per_columns[n] <- 2
# }
}

if( "unset" %in% colors_per_columns) stop('missing color assingment to some samples')

# convert colors to matrix of rgb values
col_matrix <- col2rgb(colors_per_columns, alpha = TRUE) / 255
# set alpha channel for all "primary tumor" samples to 0.1
col_matrix[4, grep("01", tissue_code)] <- 0.1

################################################################################
# create the CDF for every column and plot it
################################################################################

used_cols <- c()

pdf(pdf_out, paper="a4", width=8, height=8*1.66)

layout(matrix(c(1,2), 2, 1, byrow = TRUE),heights=c(2,4))
plot(0,type="n",xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",axes=F)
exons_gt_0 <- lengthPerExon[[1]][lengthPerExon[[1]]>0]
first_ecdf <- ecdf(exons_gt_0)
legendValues <- c(length(exons_gt_0))

plot(first_ecdf, do.points = F, verticals = F, lwd=2, main=plotMain, xlab = xlabel, ylab=ylabel, col = rgb(col_matrix[1,1], col_matrix[2,1], col_matrix[3,1], col_matrix[4,1]))
used_cols <- c(used_cols, rgb(col_matrix[1,1], col_matrix[2,1], col_matrix[3,1], col_matrix[4,1]) )
for(i in 2:dim(lengthPerExon)[2]){
  exons_gt_0 <- lengthPerExon[[i]][lengthPerExon[[i]]>0]
  current_ecdf <- ecdf(exons_gt_0)
  curr_col <- rgb(col_matrix[1,i], col_matrix[2,i], col_matrix[3,i], col_matrix[4,i])
  lines(current_ecdf, do.points = F, verticals = F, col = curr_col, lwd=2)
  used_cols <- c(used_cols, curr_col)
  legendValues <- c(legendValues, length(exons_gt_0))
}

# plot legend on other page
par(mfrow=c(1,1))
plot(0,type="n",xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",axes=F)
legend("topleft", legend = paste(colnames(lengthPerExon), " (", legendValues, ")"), col = used_cols, lwd=2, cex = 0.8)
dev.off()

