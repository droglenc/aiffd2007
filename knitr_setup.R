library(knitr)

################################################################################
## knitr custom hooks
################################################################################
knit_hooks$set(
  par1 = function(before, options, envir) {
    if (before && options$fig.show != "none") {
      par(mar=c(3.05,3.05,0.65,0.65),mgp=c(1.9,0.3,0),
          tcl=-0.2,las=1,cex.lab=0.95,cex.axis=0.9)
      if (is.list(options$par))
        do.call(par, options$par)
    }
  }
)

################################################################################
## knitr options -- figure handling
################################################################################
# set default plot options
opts_chunk$set(fig.path='Figs/',par1=TRUE)
# output look
opts_chunk$set(prompt=TRUE,comment='')


################################################################################
## r options
################################################################################
options(width=80,show.signif.stars=FALSE,continue="  ",
        str=strOptions(strict.width="cut"))
