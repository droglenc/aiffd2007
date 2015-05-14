library(FSA)          # Summarize, fitPlot, residPlot, addSigLetters, lencat, headtail
library(NCStats)      # addSigLetters
library(car)          # Anova
library(multcomp)     # glht, mcp
library(lattice)      # xyplot

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,0.5,0.5),mgp=c(2.1,0.4,0),las=1,tcl=-0.2)

d1 <- read.table("data/Box10_1.txt",header=TRUE)
str(d1)
d1$logTL <- log10(d1$TL)
d1$logWT <- log10(d1$WT)
str(d1)

d1A <- Subset(d1,POP=="A")
str(d1A)
d1B <- Subset(d1,POP=="B")
str(d1B)

lm1 <- lm(logWT~logTL,data=d1A)
anova(lm1)
summary(lm1)
confint(lm1)

fitPlot(lm1,xlab="log10(Total Length)",ylab="log10(Weight)",main="")

residPlot(lm1)

lm2 <- lm(logWT~logTL,data=d1B)
anova(lm2)
summary(lm2)
confint(lm2)

lens <- c(250,300,450)    # setup lengths
predA.logwt <- predict(lm1,data.frame(logTL=log10(lens)),interval="c")
10^(predA.logwt)          # back-transform to raw weights
predB.logwt <- predict(lm2,data.frame(logTL=log10(lens)),interval="c")
10^(predB.logwt)

nls1 <- nls(WT~A*TL^B,data=d1A,start=list(A=0.000001,B=3.0),trace=TRUE)

summary(nls1,corr=TRUE)
confint(nls1)

fitPlot(nls1,xlab="Total Length",ylab="Weight",main="")

residPlot(nls1)

predict(nls1,data.frame(TL=lens))

d1AB <- Subset(d1,POP!="C")
str(d1AB)

lm3 <- lm(logWT~logTL*POP,data=d1AB)
Anova(lm3,type="II")

lm3r <- lm(logWT~logTL+POP,data=d1AB)
Anova(lm3r,type="II")
summary(lm3r)
confint(lm3r)
fitPlot(lm3r,xlab="log10 Total Length",ylab="log10 Weight",main="",legend="topleft")

d1$WS <- as.numeric(NA)
d1$WS[d1$POP=="A"] <- 10^(-5.189+3.099*d1$logTL[d1$POP=="A"])  # Lotic
d1$WS[d1$POP!="A"] <- 10^(-5.192+3.086*d1$logTL[d1$POP!="A"])  # Lentic
d1$WR <- (d1$WT/d1$WS)*100
str(d1)

d1 <- lencat(~TL,data=d1,breaks=c(100,150,200,250,300),vname="GRP")
levels(d1$GRP)
headtail(d1)

d1A <- Subset(d1,POP=="A")
str(d1A)

lm1 <- lm(WR~GRP,data=d1A)
anova(lm1)
summary(lm1)
fitPlot(lm1,xlab="Total Length Category",ylab="Relative Weight",main="")

residPlot(lm1)

kruskal.test(WR~GRP,data=d1A)

lm2 <- lm(WR~POP,data=d1)
anova(lm2)
summary(lm2)
Summarize(WR~POP,data=d1,digits=2)
fitPlot(lm2,xlab="Population",ylab="Relative Weight",main="")

residPlot(lm2)

lm3 <- lm(FAT~WR,data=d1)
anova(lm3)
summary(lm3)
fitPlot(lm3,xlab="Relative Weight",ylab="Fat Composition",main="")

residPlot(lm3)

lm4 <- lm(FAT~POP,data=d1)
anova(lm4)
Summarize(FAT~POP,data=d1,digits=2)

residPlot(lm4)

mc1 <- glht(lm4,mcp(POP="Tukey"))
summary(mc1)
confint(mc1)
cld(mc1)

fitPlot(lm4,xlab="Population",ylab="Fat Composition",main="")
addSigLetters(lm4,lets=c("a","b","b"),pos=c(2,2,4))

d <- read.table("data/Box10_5.txt",header=TRUE)
str(d)

lm1 <- lm(BD~SL*food,data=d)
Anova(lm1,type="II")

lm1r <- lm(BD~SL+food,data=d)
Anova(lm1r,type="II")
summary(lm1r)
confint(lm1r)
fitPlot(lm1r,xlab="Standard Length (mm)",ylab="Body Depth (mm)",main="",legend="topleft")

residPlot(lm1r)

d6 <- read.table("data/Box10_6.txt",header=TRUE)
str(d6)
d6$logTL <- log10(d6$TL)
d6$logWT <- log10(d6$WT)
d6$K <- d6$WT/(d6$TL^3)*10000   # Fulton's condition factor
headtail(d6)                        # six random rows from the data frame

d6w <- Subset(d6,Lernaea=="with")
str(d6w)
d6wo <- Subset(d6,Lernaea=="without")
str(d6wo)

Summarize(TL~Lernaea,data=d6)
Summarize(WT~Lernaea,data=d6)
Summarize(K~Lernaea,data=d6)

lmw <- lm(logWT~logTL,data=d6w)
anova(lmw)
summary(lmw)
confint(lmw)
lmwo <- lm(logWT~logTL,data=d6wo)
anova(lmwo)
summary(lmwo)                                                                                                        
confint(lmwo)

# for Sexpr below
lmc <- lm(logWT~logTL*Lernaea,data=d6)

lmc <- lm(logWT~logTL*Lernaea,data=d6)
Anova(lmc,type="II")
summary(lmc)
fitPlot(lmc,xlab="log10 Total Length",ylab="log10 Weight",main="",legend="topleft")

residPlot(lmc,main="")

lmTL <- lm(TL~Lernaea,data=d6)
anova(lmTL)
summary(lmTL)
lmWT <- lm(WT~Lernaea,data=d6)
anova(lmWT)
summary(lmWT)
lmK <- lm(K~Lernaea,data=d6)
anova(lmK)
summary(lmK)

plot(WT~TL,data=d6,subset=Lernaea=="with",ylab="Wt(g)",xlab="TL (mm)",pch=16)
points(WT~TL,data=d6,subset=Lernaea=="without",pch=1)
legend("topleft",legend=c("Lernaea present","Lernaea absent"),pch=c(16,1),cex=0.8)
plot(logWT~logTL,data=d6,subset=Lernaea=="with",ylab="log10(Wt)",xlab="log10(TL)",pch=16)
points(logWT~logTL,data=d6,subset=Lernaea=="without",pch=1)
plot(K~TL,data=d6,subset=Lernaea=="with",ylab="K",xlab="TL (mm)",pch=16)
points(K~TL,data=d6,subset=Lernaea=="without",pch=1)

xyplot(WT~TL,group=Lernaea,data=d6,xlab="Total Length (mm)",ylab="Weight(g)",auto.key=TRUE)
xyplot(WT~TL,group=Lernaea,data=d6,xlab="Total Length (mm)",ylab="Weight(g)",pch=c(16,1),
       col="black",key=list(text=list(c("Lernaea present","Lernaea absent")),
                            points=list(pch=c(16,1),col="black")))
xyplot(K~TL,group=Lernaea,data=d6,xlab="Total Length (mm)",ylab="K",pch=c(16,1),
       col="black",key=list(text=list(c("Lernaea present","Lernaea absent")),
                            points=list(pch=c(16,1),col="black")))


# Script created at 2015-05-14 16:34:27
