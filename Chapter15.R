library(FSA)          #

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

d <- read.table("data/Box15_1.txt",header=TRUE)
str(d)

dt <- t(d[,-c(1:2)])        # transpose only the count portion
colnames(dt) <- d$SpecCode  # use species codes to label columns (variable names)
dt

res1 <- prcomp(dt,retx=TRUE,center=TRUE,scale.=TRUE)
summary(res1)
screeplot(res1)
biplot(res1)
pc1 <- res1$rotation[,1]  # first PCA loadings
pc2 <- res1$rotation[,2]
pc1s <- predict(res1)[,1] # first PCscore for each individual


# Script created at 2015-05-15 17:34:38
