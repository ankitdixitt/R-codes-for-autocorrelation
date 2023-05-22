directory <- "H:/Ankit MTech/FCM data for Ankit/Output_FCS-to-CSV_2023-05-16_04-30-14.393833"
# primers <- as.character(c(1:5))


carbonsource <- c("NCD", "NHD", "NLD", "NMD", "OHD", 
                  "OLD", "OMD", "PCD")
time <- as.character(c(0:9))
rep <- as.character(c(1:3))
binw <- 1 ## binwidth - we are simply rounding the base size to the closest integer
library(IDPmisc)
library(ggplot2)
library(ecodist)
library(plyr)
library(reshape)
library(aws)
library(ca)
library(rgl)
library(FactoMineR)
library(vcd)

install.packages("ggplot2")






### Program to plot ARISA data
# library(emdist)
library(ash)
library(flowCore)
directory <- "H:/Ankit MTech/FCM Data For Ankit/"
fset <- read.flowSet(path=directory, pattern="*.fcs", truncate_max_range = FALSE)
save(fset, file="fset.RData")
## This chunk crashes because of limited memory (at 8 GB you'd think it was 
## enough.. but no!)
## So I computed this part separately and will load it in the next chunk. 

disthelper <- function(ff1, ff2, col1, col2, bins) {
  matrix1 <- bin2(asinh(exprs(ff1)[,c(col1, col2)]), nbin=c(bins,bins))$nc
  matrix2 <- bin2(asinh(exprs(ff2)[,c(col1, col2)]), nbin=c(bins,bins))$nc
  
  matrix1 <- matrix1 / sum(matrix1)
  matrix2 <- matrix2 / sum(matrix2)
  
  list(matrix1, matrix2)
}

emdvalue <- function(ff1, ff2, col1, col2, bins) {
  matrices <- disthelper(ff1, ff2, col1, col2, bins)
  emd2d(matrices[[1]], matrices[[2]])
}

rmsd <- function(ff1, ff2, col1, col2, bins) {
  matrices <- disthelper(ff1, ff2, col1, col2, bins)
  sqrt(mean((matrices[[1]] - matrices[[2]])^2))
}







# <<flowcyt_distmats, cache=FALSE, echo=FALSE, eval=FALSE>>=
  ## The below is all commented because it crashes while reading and writing to the
  ## cache database.
  
load("fset.RData")
distMatrixFA <- matrix(0, nrow=length(fset), ncol=length(fset))
for(i in c(1:length(fset))) {
  for (j in c(1:length(fset))) {
    distMatrixFA[i,j] <- rmsd(fset[[i]], fset[[j]], "FSC-A", "SSC-A", 50)
  }
}
save(distMatrixFA, file = "FAdm.RData")
#
#


# ################################################################################# 
# distMatrixSA <- matrix(0, nrow=length(fset), ncol=length(fset))
# for(i in c(1:length(fset))) {
#   for (j in c(1:length(fset))) {
#     distMatrixSA[i,j] <- rmsd(fset[[i]], fset[[j]], "SSC-A", "AmCyan-A", 50)
#   }
# }
# save(distMatrixSA, file ="SAdm.RData")
# #
# #
# distMatrixFS <- matrix(0, nrow=length(fset), ncol=length(fset))
# for(i in c(1:length(fset))) {
#   for (j in c(1:length(fset))) {
#     distMatrixFS[i,j] <- rmsd(fset[[i]], fset[[j]], "FSC-A", "SSC-A", 50)
#   }
# }
# save(distMatrixFS, file = "FSdm.RData")
# 


################################################################################
## We will create a 3D structure based on the FSC-A, SSC-A and AmCyan-A data.

## Values in the FSC / SSC / AmCyan table go from 0-10000. We define a bin number of 50 along each dimension.
# 
# 
xrangemat <- t(matrix(rep(c(0,asinh(10000)), 2), ncol=2, byrow=T))
# 
bindat <- function(ff) {
  binning(asinh(exprs(ff[,c('FSC-A', 'SSC-A')])),
          y=NULL, nbins=c(50,50),
          xrange = xrangemat)
}
# 
array_rmsd <- function(d1, d2) {
  # Function to compute RMSD between two density arrays
  table1 <- bindat(d1)$table.freq / 10000
  table2 <- bindat(d2)$table.freq / 10000
  sqrt(mean((table1 - table2)^2))
}
# 
RMSD3d <- matrix(0, nrow=length(fset), ncol=length(fset))
for(i in 1:length(fset)) {
  print(i)
  for(j in (i+1):length(fset)) {
    RMSD3d[i,j] <- array_rmsd(fset[[i]], fset[[j]])
  }
}
# 
RMSD3d <- RMSD3d + t(RMSD3d)
save(RMSD3d, file='RMSD3d.RData')


# cardinale_correction <- TRUE ## Increase basepair lengths to synchronize


directory <- "H:/Ankit MTech/FCM data for Ankit/Output_FCS-to-CSV_2023-05-16_04-30-14.393833/"
filenames <- dir(path = directory, pattern = ".csv")
data <- ldply(paste0(directory, filenames), function(f) {
  df <- read.csv(f, colClasses=c("numeric"))
  if(nrow(df)>0) df$filename <- as.factor(f)
  df
})


# <<flowcyt_init, cache=TRUE, echo=FALSE, message=FALSE>>=
  ### Program to plot ARISA data
  # library(emdist)
library(ash)
library(flowCore)
directory <- "H:/Ankit MTech/FCM Data For Ankit/"
fset <- read.flowSet(path=directory, pattern="*.fcs")
save(fset, file="fset.RData")
## This chunk crashes because of limited memory (at 8 GB you'd think it was 
## enough.. but no!)
## So I computed this part separately and will load it in the next chunk. 

disthelper <- function(ff1, ff2, col1, col2, bins) {
  matrix1 <- bin2(asinh(exprs(ff1)[,c(col1, col2)]), nbin=c(bins,bins))$nc
  matrix2 <- bin2(asinh(exprs(ff2)[,c(col1, col2)]), nbin=c(bins,bins))$nc
  
  matrix1 <- matrix1 / sum(matrix1)
  matrix2 <- matrix2 / sum(matrix2)
  
  list(matrix1, matrix2)
}

emdvalue <- function(ff1, ff2, col1, col2, bins) {
  matrices <- disthelper(ff1, ff2, col1, col2, bins)
  emd2d(matrices[[1]], matrices[[2]])
}

rmsd <- function(ff1, ff2, col1, col2, bins) {
  matrices <- disthelper(ff1, ff2, col1, col2, bins)
  sqrt(mean((matrices[[1]] - matrices[[2]])^2))
}








# <<distmat_load, cache=TRUE, echo=FALSE>>=
load("FAdm.RData")
load('RMSD3d.RData')
load("filenames.RData")

filenames <- sub("BLAN", "NONE", filenames)
# matrix_names <- paste(substr(filenames, 3,6), "T",
# substr(filenames, 8,10), 
# sep="_")
matrix_names <- filenames
# dimnames(distMatrixFS) <- list(matrix_names, matrix_names)
dimnames(distMatrixFA) <- list(matrix_names, matrix_names)
# dimnames(distMatrixSA) <- list(matrix_names, matrix_names)
dimnames(RMSD3d) <- list(matrix_names, matrix_names)


rmsd <- function(ff1, ff2, col1, col2, bins) {
  matrices <- disthelper(ff1, ff2, col1, col2, bins)
  sqrt(mean((matrices[[1]] - matrices[[2]])^2))
}
















# <<hist_fset_distances, cache=TRUE, echo=FALSE, fig.width=7, fig.height=3, fig.align='center'>>=
  
  getdists <- function(mat) {
    within_dist <- numeric(0)
    between_dist <- numeric(0)
    for(i in c(1:(nrow(as.matrix(mat))/3))) {
      vec <- (i-1)*3 + c(1,2,3)
      submat_within <- as.matrix(mat)[vec, vec]
      dists_w <- as.numeric(as.dist(submat_within))
      ## the above 3 lines extract the elements of a 3x3 matrix for every triplicate
      within_dist <- c(within_dist, dists_w)
      
      if(i < 14) {
        submat_between <- as.matrix(mat)[vec, c(((i-1)*3+4): 
                                                  ncol(as.matrix(mat)))]
        dists_bet <- as.numeric(submat_between)
        between_dist <- c(between_dist, dists_bet)
      }
    }
    
    ## combine into data frame for ggplot density plots
    within <- data.frame(dist = within_dist)
    within$type <- "Within"
    
    between <- data.frame(dist = between_dist)
    between$type <- "Between"
    
    rbind(within, between)
    
  }

# ddfa <- getdists(distMatrixFS)
# ddfa$data <- "FSC-AmCyan"

ddfs <- getdists(distMatrixFA)
ddfs$data <- "FSC-SSC"

# ddsa <- getdists(distMatrixSA)
# ddsa$data <- "SSC-AmCyan"

dd3d <- getdists(RMSD3d)
dd3d$data <- "3D Voxel"

bigd <- rbind(ddfs, dd3d)
ggplot(bigd, aes(dist, fill=type)) + geom_density(alpha=0.2) +
  scale_x_log10() + facet_grid(~data)













# <<kst_flowcyt, cache=TRUE, echo=FALSE>>=
  kst_fun <- function(dat) {
    ks.test(subset(dat, type=="Within")$dist, subset(dat, type=="Between")$dist, 
            alternative="g")
  }
# kst_fun(ddfa)
kst_fun(ddfs)
# kst_fun(ddsa)
kst_fun(dd3d)
#' @
  # <<fset_crosscor, cache=TRUE, echo=FALSE>>=
  cor(data.frame(
    # as.numeric(distMatrixFS), 
    as.numeric(distMatrixFA), 
    # as.numeric(distMatrixSA), 
    as.numeric(RMSD3d)))











# <<autocor_fset, echo=FALSE, cache=TRUE, fig.height=2, fig.width=6, fig.align='center'>>=
autocor <- numeric(0)
ab <- unique(substr(filenames, 1,3))
# ab <- ab[-4]

extract_timeseries <- function(temp_source, mat, rep) {
  ## For a given carbon source get the matrix containing the time series data
  ## for that particular replicate
  day <- c(0:9)
  # vec <- c(paste("NONE_T_0_", rep, sep=""),
  #                 paste(c(carbonsource), c("T"), day, c(rep), sep="_"))
  
  vec <- c(paste(c(temp_source), rep(0:9, each = 3), "_", rep(1:3, 10), ".csv", sep = ""))
  
  # vec <- c(c(paste0("SLD", rep(1:5, each = 3), "_", rep(1:3, 5), ".csv")), 
  #          c(paste0("SMD", rep(1:5, each = 3), "_", rep(1:3, 5), ".csv")),
  #          c(paste0("SHD", rep(1:5, each = 3), "_", rep(1:3, 5), ".csv")),
  #          c(paste0("PCD", rep(1:5, each = 3), "_", rep(1:3, 5), ".csv")),
  #          c(paste0("NCD", rep(1:5, each = 3), "_", rep(1:3, 5), ".csv")))
  
  
  mat[vec,vec]
  
  
}

delay_vec <- function(distmat, delay) {
  vec <- numeric(0)
  for(temp_ab in ab) {
    for(rep in c(1:3)) {
      submat <- extract_timeseries(temp_ab, distmat, rep)
      diags <- diag(submat[c(1:(11-delay)), c((1+delay):11)])
      df <- data.frame(diags)
      df$temp_ab <- temp_ab
      df$replicate <- rep
      df$time <- c(1:length(diags)) + (delay-1)
      
      vec <- rbind(vec, df)
    }
  }
  vec
}

get_delays <- function(distmat) {
  mat <- numeric(0)
  for(delay in c(1,2,3,4,5,6,7,8,9)) {
    temp <- delay_vec(distmat, delay)
    temp$delay <- delay
    mat <- rbind(mat, temp)
  }
  mat
}

write.csv(RMSD3d, "RMSD3d.csv")

ggplot(get_delays(RMSD3d), aes(x=as.factor(delay), y=diags, fill=temp_ab)) +   
  geom_boxplot() + scale_y_log10(name="Dissimilarity (RMSD 3D)") + facet_grid(~temp_ab) + 
  scale_x_discrete(name="Delay (d)") +
  scale_y_log10()

