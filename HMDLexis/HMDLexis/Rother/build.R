if (system("hostname",intern=TRUE) %in% c("triffe-N80Vm", "tim-ThinkPad-L440")){
	# if I'm on the laptop
	wd <- "/home/tim/git/HMDLexis/HMDLexis/HMDLexis"
} else {
	# in that case I'm on Berkeley system, and other people in the dept can run this too
	wd <- paste0("/data/commons/",system("whoami",intern=TRUE),"/git/HMDLexis/HMDLexis/HMDLexis")
}

library(devtools)
document(wd)

#install.packages("pspline")
#install.packages("pracma")
install_github("timriffe/TimUtils/TimUtils")
library(TimUtils)
parent.path <- "/data/commons/triffe/git/HMDLexis/HMDLexis"
IncrementVersion(file.path(parent.path ,"HMDLexis"),"5","2012-12-05")

