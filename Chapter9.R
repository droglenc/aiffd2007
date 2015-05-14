library(FSA)           # Summarize, fitPlot, addSigLetters, residPlot
library(NCStats)       # chisqPostHoc
library(lattice)       # bwplot, xyplot
library(multcomp)      # glht, mcp
library(nlme)          # lme
library(pgirmess)      # kruskalmc
library(TeachingDemos) # chisq.detail

setwd("C:/aaaWork/Web/fishR/BookVignettes/aiffd2007")

options(width=90,continue=" ",show.signif.stars=FALSE,
        contrasts=c("contr.sum","contr.poly"))
par(mar=c(3.5,3.5,1,1),mgp=c(2.1,0.4,0),tcl=-0.2)

d1 <- read.table("data/Box9_1.txt",header=TRUE)
str(d1)

lm1 <- lm(length~lake,data=d1)
anova(lm1)

mc1 <- glht(lm1,mcp(lake="Tukey"))
summary(mc1)
confint(mc1)
cld(mc1)

fitPlot(lm1,ylab="Total Length (mm)",xlab="Lake",main="")
addSigLetters(lm1,lets=c("c","a","b"),pos=c(2,2,4))

Summarize(length~lake,data=d1,digits=2)

d1I <- Subset(d1,lake=="Island")     # only Island
d1M <- Subset(d1,lake=="Mitchell")   # only Mitchell
d1T <- Subset(d1,lake=="Thompson")   # only Thompson

ks.test(d1M$length,d1T$length)
ks.test(d1I$length,d1T$length)
ks.test(d1I$length,d1M$length)

hist(length~lake,data=d1,xlab="Total Length (mm)")

plot(ecdf(d1I$length),xlim=c(100,240),verticals=TRUE,pch=".",main="",
     xlab="Total Length (mm)",lwd=2)
plot(ecdf(d1M$length),col="blue",verticals=TRUE,pch=".",lwd=2,add=TRUE)
plot(ecdf(d1T$length),col="red",verticals=TRUE,pch=".",lwd=2,add=TRUE)
legend(100,0.99,legend=c("Island","Mitchell","Thopmson"),col=c("black","blue","red"),
       pch=".",lwd=2,lty=1,cex=0.75)

kruskal.test(d1$length,d1$lake)

kruskalmc(d1$length,d1$lake)

d5 <- matrix(c(85,77,44,124,34,251),nrow=3,byrow=TRUE)
colnames(d5) <- c("Q","S-Q")
rownames(d5) <- c("1996","1997","1998")
d5

( chi1 <- chisq.test(d5) )
chi1$expected
chi1$residuals

chisq.detail(d5)

chisqPostHoc(chi1,digits=6)

d6 <- read.table("data/Box9_6.txt",header=TRUE)
str(d6)
d6$LOGIT <- log((d6$PREF+0.5)/(d6$QUAL+0.5))
d6

Summarize(LOGIT~METHOD,data=d6,digits=4)

lm1 <- lm(LOGIT~METHOD,data=d6,weights=QUAL)
anova(lm1)

d7 <- read.table("data/Box9_7.txt",header=TRUE)
str(d7)
levels(d7$sizegrp)
d7$fsite <- factor(d7$site)
d7$fyear <- factor(d7$year)

d7$period <- "APRE"                  # initially fill completely with "APRE"
d7$period[d7$year>1990] <- "BPOST"   # then replace post-1990 with "BPOST"
d7$period <- factor(d7$period)       # explicitly make a factor
view(d7)
d7a <- Subset(d7,sizegrp!="age0")

d7b <- reshape(d7a,direction="wide",timevar="sizegrp",idvar=c("site","year","fyear","period"))
view(d7b)

d7b$undert <- d7b$count.und+0.5
d7b$slott <- d7b$count.slot+0.5
d7b$total <- d7b$undert + d7b$slott
d7b$LOGIT <- log(d7b$slott/d7b$undert)
head(d7b)      # first 6 rows

shapiro.test(d7b$LOGIT[d7b$period=="APRE"])
shapiro.test(d7b$LOGIT[d7b$period=="BPOST"])

qqnorm(d7b$LOGIT[d7b$period=="APRE"],main="Pre-Regulation")
qqnorm(d7b$LOGIT[d7b$period=="BPOST"],main="Post-Regulation")

me1 <- lme(LOGIT~period,random=~1|fyear,weights=~total,correlation=corAR1(0,~fsite),data=d7b)

library(lattice)
dat <- read.table("data/http://www.ncfaculty.net/dogle/me_ex.txt", head=TRUE)
names(dat) <- c('yr','site','per','tot','logit')
dat <- transform(dat, yrf=factor(yr), site=factor(site))
d0 <- dat

#d0 <- groupedData(logit~per|yr, d0)
bwplot(logit~per|site, d0)
xyplot(logit~yr|site, d0, group=per, auto.key=TRUE)
m2 <- lme(logit~per, data=d0, random=~1|site, weights=~tot, correlation=corAR1(form=~yr|site))

m2a <- lme(logit~0+per, data=d0, random=~1|site, weights=~tot, correlation=corAR1(form=~yr|site))

m3a <- lme(logit~0+per, data=d0, random=~1|yr, weights=~tot, correlation=corAR1(form=~yr|site))


# Script created at 2015-05-14 17:45:03
