library(FSA)          # Subset, Summarize, residPlot
library(car)          # Anova, recode
library(plyr)         # ddply

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

d2 <- read.table("data/Box12_2.txt",header=TRUE)
str(d2)

d2$tempA <- NA
d2$tempA[d2$temp>6 & d2$temp<9] <- 6.9
d2$tempA[d2$temp>21 & d2$temp<24] <- 22.4
d2$tempA[d2$temp>28 & d2$temp<31] <- 29.6
d2$tempA <- factor(d2$tempA)
str(d2)

d2$wtA <- NA
d2$wtA[d2$wt<77] <- 38
d2$wtA[d2$wt>130 & d2$wt<730] <- 403
d2$wtA[d2$wt>800] <- 1567
d2$wtA <- factor(d2$wtA)
str(d2)

d2$logwt <- log10(d2$wt)
d2$logCmax <- log10(d2$Cmax)
lm.low <- lm(Cmax~wt,data=d2,subset=tempA=="6.9")
coef(lm.low)
lm.med <- lm(logCmax~logwt,data=d2,subset=tempA=="22.4")
coef(lm.med)
lm.hi  <- lm(logCmax~logwt,data=d2,subset=tempA=="29.6")
coef(lm.hi)

# schematic plot
plot(Cmax~wt,data=d2,type="n",xlab="Wet weight (g)",ylab="Cmax(g/g/d)")
# add points specific to each group
points(Cmax~wt,data=d2,subset=tempA=="6.9",pch=19,col="red")
points(Cmax~wt,data=d2,subset=tempA=="22.4",pch=19,col="blue")
points(Cmax~wt,data=d2,subset=tempA=="29.6",pch=19,col="black")
# add the fitted curve specific to each group
curve(coef(lm.low)[1]+coef(lm.low)[2]*x,min(d2$wt[d2$tempA=="6.9"],na.rm=TRUE),
      max(d2$wt[d2$tempA=="6.9"],na.rm=TRUE),add=TRUE,col="red",lwd=2)
curve(10^(coef(lm.med)[1])*x^coef(lm.med)[2],min(d2$wt[d2$tempA=="22.4"],na.rm=TRUE),
      max(d2$wt[d2$tempA=="22.4"],na.rm=TRUE),add=TRUE,col="blue",lwd=2)
curve(10^(coef(lm.hi)[1])*x^coef(lm.hi)[2],min(d2$wt[d2$tempA=="29.6"],na.rm=TRUE),
      max(d2$wt[d2$tempA=="29.6"],na.rm=TRUE),add=TRUE,col="black",lwd=2)
# add a legend
legend("topright",legend=c("6.9C","22.4C","29.6C"),pch=19,col=c("red","blue","black"),
       lty=1,lwd=2)

lm.sm <- lm(Cmax~temp+I(temp^2),data=d2,subset=wtA=="38")
coef(lm.sm)
lm.int <- lm(Cmax~temp+I(temp^2),data=d2,subset=wtA=="403")
coef(lm.int)
lm.lrg  <- lm(Cmax~temp+I(temp^2),data=d2,subset=wtA=="1567")
coef(lm.lrg)

# The plot
plot(Cmax~temp,data=d2,type="n",xlab="Temperature (C)",ylab="Cmax(g/g/d)")
points(Cmax~temp,data=d2,subset=wtA=="38",pch=19,col="red")
points(Cmax~temp,data=d2,subset=wtA=="403",pch=19,col="blue")
points(Cmax~temp,data=d2,subset=wtA=="1567",pch=19,col="black")
curve(coef(lm.sm)[1]+coef(lm.sm)[2]*x+coef(lm.sm)[3]*x^2,min(d2$temp,na.rm=TRUE),
      max(d2$temp,na.rm=TRUE),add=TRUE,col="red",lwd=2)
curve(coef(lm.int)[1]+coef(lm.int)[2]*x+coef(lm.int)[3]*x^2,min(d2$temp,na.rm=TRUE),
      max(d2$temp,na.rm=TRUE),add=TRUE,col="blue",lwd=2)
curve(coef(lm.lrg)[1]+coef(lm.lrg)[2]*x+coef(lm.lrg)[3]*x^2,min(d2$temp,na.rm=TRUE),
      max(d2$temp,na.rm=TRUE),add=TRUE,col="black",lwd=2)
legend("topleft",legend=c("38 g","403 g","1567 g"),pch=19,col=c("red","blue","black"),
       lty=1,lwd=2)

d2a <- Subset(d2,Cmax>0)

lm1 <- lm(logCmax~logwt*temp,data=d2a)
anova(lm1)
Anova(lm1,type="III")

lm2 <- lm(logCmax~logwt+temp+I(temp^2),data=d2a)
anova(lm2)
Anova(lm2,type="III")
summary(lm2)

d2a$pCmax <- 10^lm2$fitted
d2a$resids2 <- d2a$Cmax-d2a$pCmax
plot(resids2~temp,data=d2a,type="n",xlab="Temperature",ylab="Residuals")
points(resids2~temp,data=d2a,subset=wtA=="38",col="red",pch=19)
points(resids2~temp,data=d2a,subset=wtA=="403",col="blue",pch=19)
points(resids2~temp,data=d2a,subset=wtA=="1567",col="black",pch=19)
abline(v=c(5,10,15,20,24,28),col="gray90",lty=3)
legend("topleft",legend=c("38 g","403 g","1567 g"),pch=19,col=c("red","blue","black"),
       lty=1,lwd=2)

d2a$tempB <- recode(d2a$temp,"0:5='3'; 6:10='8'; 11:15='13'; 16:20='18'; 
                    21:23.5='22'; 24:28='26'; else='29'") 
sumcv <- ddply(d2a,c("wtA","tempB"),function(df) sd(df$pCmax)/mean(df$pCmax)*100)
names(sumcv)[3] <- "CV"
str(sumcv)

plot(CV~tempB,data=sumcv,type="n",xlab="Temperature (C)",ylab="CV")
points(CV~tempB,data=sumcv,subset=wtA==38,type="b",pch=19,col="red",lwd=2)
points(CV~tempB,data=sumcv,subset=wtA==403,type="b",pch=19,col="blue",lwd=2)
points(CV~tempB,data=sumcv,subset=wtA==1567,type="b",pch=19,col="black",lwd=2)
legend("topright",legend=c("38 g","403 g","1567 g"),pch=19,col=c("red","blue","black"),
       lty=1,lwd=2)

plot(lm2$residuals~lm2$fitted,pch=19,xlab="Fitted Values",ylab="Residuals")


# Script created at 2015-05-15 18:10:28
