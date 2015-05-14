library(FSA)          # Subset, fitPlot, ks2d2, ks2d2p, headtail
library(car)          # Anova
library(ggplot2)      # qplot, facet_grid
library(Hmisc)        # rcorr
library(lattice)      # xyplot
library(MASS)         # manova
library(quantreg)     # rq et al.

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

d1 <- matrix(c(18,9,18,25),nrow=2,byrow=TRUE)
colnames(d1) <- c("J.present","J.absent")
rownames(d1) <- c("I.present","I.absent")
d1
addmargins(d1)                                 # for display only

( g1 <- g.test(d1) )

( chi1 <- chisq.test(d1,correct=FALSE) )

(d1[1,1]*d1[2,2])/(d1[1,2]*d1[2,1])

d2 <- read.table("data/Box11_2.txt",header=TRUE)
str(d2)

xyplot(abundance~freqocc|pred1,data=d2,pch=19,layout=c(1,4),as.table=TRUE,
       xlab="Frequency of Occurrence",ylab="Prey-Specific Abundance")

qplot(freqocc,abundance,data=d2,xlab="Frequency of Occurrence",
      ylab="Prey-Specific Abundance") + facet_grid(pred1~.)

d3 <- read.table("data/Box11_3.txt",header=TRUE)
str(d3)

plot(PREY.len~LMB.len,data=d3,pch=19,xlab="Largemouth Bass Length (mm)",
     ylab="Prey Fish Length (mm)")

# for Sexpr below
qr1 <- rq(PREY.len~LMB.len,data=d3,tau=0.5)

qr1 <- rq(PREY.len~LMB.len,data=d3,tau=0.5)
summary(qr1,covariance=TRUE)
summary(qr1,covariance=TRUE,se="boot")

taus <- c(0.05,0.25,0.50,0.75,0.95)
qr2 <- rq(PREY.len~LMB.len,data=d3,tau=taus)
summary(qr2)
coef(qr2)

plot(PREY.len~LMB.len,data=d3,pch=19,xlab="Largemouth Bass Length (mm)",
     ylab="Prey Fish Length (mm)")
for (i in 1:length(taus)) abline(coef=coef(qr2)[,i],col=gray(taus[i]),lwd=2)

taus <- (5:95)/100
qr3 <- rq(PREY.len~LMB.len,data=d3,tau=taus)
sqr3 <- summary(qr3)
plot(qr3)

plot(sqr3)

d5 <- read.table("data/Box11_5.txt",header=TRUE)
str(d5)

#Split the data frame into the two scenarios
d5A <- Subset(d5,scenario=="A")
d5B <- Subset(d5,scenario=="B")

#Split each scenario into the fish present and vegetation present distributions,
d5Af <- Subset(d5A,site.char=="fish.present")
d5Av <- Subset(d5A,site.char=="vegetation.present")
d5Bf <- Subset(d5B,site.char=="fish.present")
d5Bv <- Subset(d5B,site.char=="vegetation.present")

plot(coord2~coord1,data=d5Av,pch=2,cex=1.5,xlab="Coordinate 1",ylab="Coordinate 2"
     ,xlim=c(1,20),ylim=c(1,4))
points(coord2~coord1,data=d5Af,pch=19,cex=1)
plot(coord2~coord1,data=d5Bv,pch=2,cex=1.5,xlab="Coordinate 1",ylab="Coordinate 2",
     xlim=c(1,20),ylim=c(1,4))
points(coord2~coord1,data=d5Bf,pch=19,cex=1)

( scenA <- ks2d2(d5Av$coord1,d5Av$coord2,d5Af$coord1,d5Af$coord2) )
plot(scenA)

( scenAp <- ks2d2p(scenA,B=1000) )
plot(scenAp)

( scenB <- ks2d2(d5Bv$coord1,d5Bv$coord2,d5Bf$coord1,d5Bf$coord2) )
plot(scenB)

( scenBp <- ks2d2p(scenB,B=1000) )
plot(scenBp)

d6 <- read.table("data/Box11_6.txt",header=TRUE)
str(d6)

# base empty figure
plot(diet~length,data=d6,type="n",xlab="Predator Length (mm)", ylab="Stomach Contents Weight (mg)")
# now add points
points(diet~length,data=d6,subset=time=="evening",pch=19,col="black")
points(diet~length,data=d6,subset=time=="morning",pch=1,col="blue")
points(diet~length,data=d6,subset=time=="noon",pch=3,col="red")
legend("topleft",legend=c("evening","morning","noon"),pch=c(19,1,3),col=c("black","blue","red"))

# for Sexpr below
lm1 <- lm(diet~length*time,data=d6)
lm2 <- lm(diet~length+time,data=d6)

lm1 <- lm(diet~length*time,data=d6)
Anova(lm1,type="III")

lm2 <- lm(diet~length+time,data=d6)
Anova(lm2,type="III")

fitPlot(lm2,xlab="Predator Length (mm)",ylab="Stomach Contents Weight (mg)",legend="topleft")

d7 <- read.table("data/Box11_7.txt",header=TRUE)
str(d7)

Y1 <- as.matrix(d7[,c("chiro","amph","odon")])

man1 <- manova(Y1~d7$lake)
summary(man1,test="Wilks")
summary(man1,test="Pillai")
summary(man1,test="Hotelling-Lawley")
summary(man1,test="Roy")

summary.aov(man1)

Y2 <- as.matrix(d7[,c("amph","odon","zoo")]) # put zoo in to get last indiv ANOVAs
man2 <- manova(Y2~d7$lake)
summary(man2,test="Wilks")                   # show that the results are the same
summary.aov(man2)

( n <- nrow(Y1) )
( E1 <- (n-1)*cov(man1$residuals) )     # Error SSCP; diagonals are indiv SSerror
( H1 <- (n-1)*cov(man1$fitted.values) ) # Treatment SSCP; diagonals are indiv SStreats

d8 <- read.table("data/Box11_8.txt",header=TRUE)
d8

glm1 <- glm(number~prey*stage*habitat,data=d8,family=poisson)
summary(glm1)
anova(glm1,test="Chisq")
Anova(glm1,type="III")

glm2 <- glm(number~prey+stage+habitat+prey:stage+prey:habitat+stage:habitat,data=d,family=poisson)
Anova(glm2,type="III")

glm3 <- glm(number~prey+stage+habitat+prey:habitat,data=d,family=poisson)
Anova(glm3,type="III")

glm3a <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey!="amphipods",family=poisson)
( aov3a <- Anova(glm3a,type="III") )

glm3c <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey!="chironomids",family=poisson)
( aov3c <- Anova(glm3c,type="III") )

glm3m <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey!="mayfiles",family=poisson)
( aov3m <- Anova(glm3m,type="III") )

glm3o <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey!="ostracods",family=poisson)
( aov3o <- Anova(glm3o,type="III") )

data.frame(prey=levels(d$prey),pvalue.int=round(c(aov3a["prey:habitat","Pr(>Chisq)"],aov3c["prey:habitat","Pr(>Chisq)"],aov3m["prey:habitat","Pr(>Chisq)"],aov3o["prey:habitat","Pr(>Chisq)"]),4))

glm3ac <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey %nin% c("amphipods","chironomids"),family=poisson)
aov3ac <- Anova(glm3ac,type="III")
glm3am <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey %nin% c("amphipods","mayflies"),family=poisson)
aov3am <- Anova(glm3am,type="III")
glm3ao <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey %nin% c("amphipods","ostacods"),family=poisson)
aov3ao <- Anova(glm3ao,type="III")
glm3cm <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey %nin% c("chironomids","mayflies"),family=poisson)
aov3cm <- Anova(glm3cm,type="III")
glm3co <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey %nin% c("chironomids","ostracods"),family=poisson)
aov3co <- Anova(glm3co,type="III")
glm3mo <- glm(number~prey+stage+habitat+prey:habitat,data=d,subset=prey %nin% c("mayflies","ostracods"),family=poisson)
aov3mo <- Anova(glm3mo,type="III")

data.frame(prey=c("amph,chiro","amph,may","amph,ostra","chiro,may","chiro,ostra","may,ostra"),pvalue.int=round(c(aov3ac["prey:habitat","Pr(>Chisq)"],aov3am["prey:habitat","Pr(>Chisq)"],aov3ao["prey:habitat","Pr(>Chisq)"],aov3cm["prey:habitat","Pr(>Chisq)"],aov3co["prey:habitat","Pr(>Chisq)"],aov3mo["prey:habitat","Pr(>Chisq)"]),4))

d9 <- read.table("data/Box11_9.txt",header=TRUE)
str(d9)
# add proportions
d9$pAmph <- d9$amph/d9$ttlwt
d9$pLFish <- d9$lfish/d9$ttlwt
d9$pDipt <- d9$dipt/d9$ttlwt
d9$pMayfly <- d9$mayfly/d9$ttlwt
# add logs of proportions
amt <- 0.00001                    # small amount to adjust for zeroes
d9$lpAmph <- log(d9$pAmph+amt)
d9$lpLFish <- log(d9$pLFish+amt)
d9$lpDipt <- log(d9$pDipt+amt)
d9$lpMayfly <- log(d9$pMayfly+amt)
# what does it all look like now
str(d9)
headtail(d9)

d9a <- d9[,c("bluegill","lpAmph","lpLFish","lpDipt","lpMayfly")]
str(d9a)

# for Sexpr below
pc1 <- princomp(d9a[,-1])
res <- summary(pc1)$sd^2
varexpl <- cumsum(res)/sum(res)

pc1 <- princomp(d9a[,-1])                    # ignore first column of bluegill IDs
summary(pc1)
screeplot(pc1,type="lines")

loadings(pc1)
biplot(pc1,xlabs=d9a$bluegill)

rcorr(cbind(d9$wt,pc1$scores[,1:2]))
plot(pc1$scores[,1]~d9$wt,pch=19,xlab="Bluegill Weight",ylab="PC1")

d10 <- read.table("data/Box11_10.txt",header=TRUE)
str(d10)
d10$P <- d10$prey.final/d10$prey.init
d10

( sums <- tapply(d10$P,d10$trial,sum) )
( sums <- rep(sums,each=2))
d10$alpha <- d10$P/sums
d10

Summarize(alpha~tag.color,data=d10,digits=2)
leveneTest(alpha~tag.color,data=d10)            # test for equal variances first
t.test(alpha~tag.color,data=d10,var.equal=TRUE)

d11 <- read.table("data/Box11_11.txt",header=TRUE)
str(d11)

plot(tp~len,data=d11,pch=19,xlab="Walleye Length (mm)",ylab="Trophic Position")

lm1 <- lm(tp~len,data=d11)
summary(lm1)
fitPlot(lm1,xlab="Walleye Length (mm)",ylab="Trophic Position",main="")

resids <- lm1$residuals
mean(abs(resids))

fitPlot(lm1,xlab="Walleye Length (mm)",ylab="Trophic Position",main="")
abline(coef=c(2.797,0.001445),lwd=2,lty=3,col="red")

resids2 <- d11$tp-(2.797+0.001445*d11$len)
mean(abs(resids2))


# Script created at 2015-05-14 18:25:14
