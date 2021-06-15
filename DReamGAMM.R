###########################################
# Set up
###########################################

library(data.table)
library(dplyr)

setwd("D:/Plyush/rats/for SA/reduced/")

###########################################
# Prepare data set
###########################################

ds_pos <- as.data.frame(fread(input = "xcms_after_IPO_MVI_XGBOOST_filter_by_QC_pos.csv", header=T))
rownames(ds_pos) <- ds_pos$V1
pos_md <- ds_pos[,c(1:6)]
ds_pos <- ds_pos[,-c(1:6)]
ds_pos <- as.data.frame(sapply(ds_pos, as.numeric))
ds_neg <- as.data.frame(fread(input = "xcms_after_IPO_MVI_XGBOOST_filter_by_QC_neg.csv", header=T))
rownames(ds_neg) <- ds_neg$V1
neg_md <- ds_neg[,c(1:6)]
ds_neg <- ds_neg[,-c(1:6)]
ds_neg <- as.data.frame(sapply(ds_neg, as.numeric))

identical(pos_md, neg_md) # check identical meta-data

ds_comb <- as.data.frame(cbind(ds_pos, ds_neg)) # combine features from both polarities
ds_comb_md <- as.data.frame(cbind(neg_md, ds_comb)) # combine all features with metadata

# filter by group
idx_del <- which(ds_comb_md$group == "QC")
ds_comb_filt <- ds_comb[-idx_del,]
ds_comb_md_filt <- ds_comb_md[-idx_del,]

# filter by day
idx_del1 <- which(ds_comb_md_filt$day == 27)
idx_del2 <- which(ds_comb_md_filt$day == 28)
idx_del <- c(idx_del1, idx_del2)
ds_comb_filt <- ds_comb_filt[-idx_del,]
ds_comb_md_filt <- ds_comb_md_filt[-idx_del,]

###########################################
# GAMM
###########################################

library(gamm4)
library(pbapply)

# data
dat <- ds_comb_md_filt
colnames(dat)[4] <- "rat" # adjust to your data
n_meta <- 6 # adjust to your data
colnames(dat)[-c(1:n_meta)] <- paste0("X", c(1:ncol(dat))) # adjust to your data
dat$rat <- as.factor(dat$rat)
dat$group <- as.factor(dat$group)
dat$day <- as.numeric(dat$day)

# perform
n_start <- 7 # adjust to your data
dss <- lapply(n_start:ncol(dat), function(y) as.data.frame(dat[, c(y, 4, 5, 6)])) # adjust to your data
dss <- lapply(1:length(dss), function(y) {colnames(dss[[y]])[1] <-"Y" 
                                          return(dss[[y]])}) # adjust to your data
pboptions(style = 3, char = ">")
gamm_fit <- pblapply(1:length(dss), function(x) gamm4(Y~s(day)+group, random=~(1|rat), data = dss[[x]])) # adjust to your data
gamm_res <- pblapply(1:length(gamm_fit), function(x) summary(gamm_fit[[x]]$gam))
gamm_pval <- as.data.frame(sapply(1:length(gamm_res), function(x) p.adjust(gamm_res[[x]][["s.table"]][4], method = "BH"))) # adjust to your data
rownames(gamm_pval) <- colnames(dat)[-c(1:n_meta)]
g_t <- 0.0
gamm_ind <- which(gamm_pval == g_t)
gamm <- colnames(ds_comb_filt)[gamm_ind]
ds_gamm <- ds_comb_filt[,gamm]
ds_gamm_f <- as.data.frame(cbind(Label = as.factor(dat$day), ds_gamm))

save(gamm_fit, file = "gamm_fit xgb.RData") # save modeling results

###########################################
# Fold-Change
###########################################

FOLD.CHANGE.MG <- function(x, f, aggr_FUN = colMeans, combi_FUN = {function(x,y) "/"(x,y)}){
  x <- log2(x)
  f <- as.factor(f)
  i <- split(1:nrow(x), f)
  x <- sapply(i, function(i){ aggr_FUN(x[i,])})
  x <- t(x)
  j <- combn(levels(f), 2)
  ret <- combi_FUN(x[j[1,],], x[j[2,],])
  rownames(ret) <- paste(j[1,], j[2,], sep = '/')
  t(ret)
}

fdr <- FOLD.CHANGE.MG(ds_gamm_f[,-1], ds_gamm_f[,1])
sel_grs <- paste0("0/", levels(ds_gamm_f$Label)[-1])
fdr_sel <- fdr[, sel_grs]
fdr_mean <- apply(abs(fdr_sel),1, mean, na.rm=T)
thr <- 1.0
fdr_thr <- names(which(fdr_mean > thr)) 
ds_gamm_fc <- as.data.frame(cbind(Label = ds_gamm_f$Label, ds_gamm_f[,fdr_thr]))

###########################################
# Dose-Response
###########################################

library(DRomics)
library(NormalizeMets)

######################### log data
ds_gamm_fc1 <- ds_gamm_fc
ds_gamm_fc1[,-1] <- LogTransform(ds_gamm_fc[,-1])$featuredata
ds_t <- as.data.frame(t(ds_gamm_fc1))
ds_t <- as.data.frame(sapply(ds_t, as.numeric))
ds_t <- as.data.frame(cbind(colnames(ds_gamm_fc), ds_t))
o <- continuousomicdata(ds_t)
s_quad <- itemselect(o, select.method = "quadratic", FDR = 0.05)
f <- drcfit(s_quad, progressbar = T)
res_dr <- f$fitres
plot(f) 

######################### raw data
ds_t1 <- as.data.frame(t(ds_gamm_fc))
ds_t1 <- as.data.frame(sapply(ds_t1, as.numeric))
ds_t1 <- as.data.frame(cbind(colnames(ds_gamm_fc), ds_t1))
o1 <- continuousomicdata(ds_t1)
s_quad1 <- itemselect(o1, select.method = "quadratic", FDR = 0.05)
f1 <- drcfit(s_quad1, progressbar = T)
res_dr1 <- f1$fitres
plot(f1)

# create final dataset
ds_gamm_fc_dr <- as.data.frame(cbind(Label = ds_gamm_fc$Label, ds_gamm_fc[, colnames(ds_gamm_fc[-1])[f[["fitres"]]$irow]])) # or f
met_id <- sapply(1:ncol(ds_gamm_fc_dr[,-1]), function(x) which(colnames(ds_comb_filt)==colnames(ds_gamm_fc_dr[,-1])[x]))
fwrite(ds_gamm_fc_dr, "ds_gamm_fc_dr.csv", row.names = T)

###########################################
# Box plots
###########################################

library(reshape2)
library(ggplot2)

df.m <- melt(ds_gamm_fc_dr, id.var = "Label") # reshape data frame

p <- ggplot(data = df.m, aes(x=variable, y=value)) + xlab("") + ylab("") +
  geom_boxplot(aes(fill=Label)) + theme(legend.position="bottom") + theme_classic() + scale_x_discrete(labels=c("")) 

(pp <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="bottom")) 

###########################################
# Scatter plots
###########################################

df.m <- melt(ds_gamm_fc_dr, id.var = "Label")
df.m <- cbind(1:nrow(df.m), df.m)
colnames(df.m)[1] <- "Patient"
p <- ggplot(data = df.m, aes(x=Label, y=value)) + xlab("") + ylab("") +
  geom_point(aes(colour=Label)) + theme(legend.position="bottom") + theme_classic()

(b <- p + facet_wrap( ~ variable, scales="free") + theme_classic() + theme(legend.position="none")) 

###########################################
# PCA
###########################################

library(FactoMineR)
library(factoextra)

pca.ds <- PCA(ds_gamm_fc_dr[,-1], scale.unit = T, graph = F)
fviz_pca_ind(pca.ds, col.ind=as.numeric(ds_gamm_fc_dr[,1]), geom = "point", 
             gradient.cols = c("red", "lightblue", "darkblue" ), legend.title = "Time", title = "")

###########################################
# End of the script
###########################################