suppressPackageStartupMessages(library("optparse"))

################################################################################
# define and load arguments/params
################################################################################

option_list <- list(
  make_option(c("-f", "--file"), type="character", help="name of the input file"),
  make_option("--main", type="character", default="Distribution avg exon length", help="title of the plot [default: %default]"),
  make_option("--xlab", type="character", default="average length [% of maximum]", help="x-axis label of the plot [default: %default]"),
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
  opt$file = ""
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
################################################################################

input <- opt$file
plotMain <- opt$main
xlabel <- opt$xlab
ylabel <- opt$ylab
pdf_out <- opt$pdf

# load input table
lengthPerExon <- read.table(input, head=T, row.names = 1, check.names = FALSE)

################################################################################
# define the set of used colors and the color for each column
################################################################################

# the number of needed colors equals the number of distinct patients
nr_colors <- length(colnames(lengthPerExon))

# use johnbaums' code of iwanthue to get distinct colors
colors <- iwanthue(nr_colors)

# preset the color per column as unset
colors_per_columns <- rep("unset", dim(lengthPerExon)[2])
# lty_per_columns <- rep(2,dim(lengthPerExon)[2])

sample_ids <- colnames(lengthPerExon)

# define the color for each column (based on the sample id)
tmp_cnt <- 0
samples_visited_already <- c()
for (n in 1:dim(lengthPerExon)[2]) {
if( sample_ids[n] %in% samples_visited_already ) {
    colors_per_columns[n] <- colors[ which(samples_visited_already == sample_ids[n]) ]
} else {
  tmp_cnt <- tmp_cnt + 1
  colors_per_columns[n] <- colors[tmp_cnt]
  samples_visited_already[tmp_cnt] <- sample_ids[n]
}
}

if( "unset" %in% colors_per_columns) stop('missing color assingment to some samples')

# convert colors to matrix of rgb values
col_matrix <- col2rgb(colors_per_columns, alpha = TRUE) / 255

################################################################################
# create the CDF for every column and plot it
################################################################################

used_cols <- c()

pdf(pdf_out, paper="a4", width=8, height=8*1.66)

layout(matrix(c(1,2), 2, 1, byrow = TRUE),heights=c(2,4))
plot(0,type="n",xlim=c(0,100),ylim=c(0,100),xlab="",ylab="",axes=F)

# Enable empty plot
if( length(lengthPerExon[[1]][lengthPerExon[[1]]>0]) == 0 ) {
    exons_gt_0 <- rep(0,length(colnames(lengthPerExon)))
} else {
    exons_gt_0 <- lengthPerExon[[1]][lengthPerExon[[1]]>0]
}

first_ecdf <- ecdf(exons_gt_0)
legendValues <- c(length(exons_gt_0))

plot(first_ecdf, do.points = F, verticals = F, lwd=2, main=plotMain, xlab = xlabel, ylab=ylabel, col = rgb(col_matrix[1,1], col_matrix[2,1], col_matrix[3,1], col_matrix[4,1]))
used_cols <- c(used_cols, rgb(col_matrix[1,1], col_matrix[2,1], col_matrix[3,1], col_matrix[4,1]) )
for(i in 2:dim(lengthPerExon)[2]){
    
    # Enable empty plot
    if( length(lengthPerExon[[i]][lengthPerExon[[i]]>0]) == 0 ) {
        exons_gt_0 <- rep(0,length(colnames(lengthPerExon)))
    } else {
        exons_gt_0 <- lengthPerExon[[i]][lengthPerExon[[i]]>0]
    }
  current_ecdf <- ecdf(exons_gt_0)
  curr_col <- rgb(col_matrix[1,i], col_matrix[2,i], col_matrix[3,i], col_matrix[4,i])
  lines(current_ecdf, do.points = F, verticals = F, col = curr_col, lwd=2)
  used_cols <- c(used_cols, curr_col)
  legendValues <- c(legendValues, length(exons_gt_0))
}

# plot legend on top left
legend("topleft", legend = paste(colnames(lengthPerExon), " (", legendValues, ")"), col = used_cols, lwd=2, cex = 0.8, bg = "white")
dev.off()

