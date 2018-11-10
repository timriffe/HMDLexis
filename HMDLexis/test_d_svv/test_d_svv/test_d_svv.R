## test routine d_svv.R in LexisLDB package
HMDLexis.loc <-  file.path("/data/commons/boe/HMDLexis.git/HMDLexis/HMDLexis")
library(devtools)
load_all(HMDLexis.loc)

setwd("/data/commons/boe/HMDLexis.git/HMDLexis/test:d_svv")

## run Matlab on the associated ldb_bgr.m script to get pre/post deaths
## d_svv.m  only called for AUS, BEL (Old), AUS, FRACNP, ITA, NLD

CTY<- 'NLD'
library(R.matlab)

load_all('/data/commons/triffe/git/LDButils/LDButils/LDButils')
deaths.pre <- loadMdeaths2R(paste0('./',CTY,'deaths_pre.mat'), PopName = CTY)
deaths.pre <- resortDeaths(deaths.pre)

deaths.post <- loadMdeaths2R(paste0('./',CTY,'deaths_post.mat'), PopName = CTY)
deaths.post <- resortDeaths(deaths.post)

deaths.rpost <- d_svv(deaths.pre, reproduce.matlab = TRUE)
deaths.rpost <- resortDeaths(deaths.rpost)

write.csv(deaths.pre, paste0(CTY,'deaths.pre.csv'), row.names=FALSE)
write.csv(deaths.post, paste0(CTY,'deaths.post.csv'), row.names=FALSE)
write.csv(deaths.rpost,paste0(CTY,'deaths.rpost.csv'), row.names=FALSE)

## refactoring : util function to deflate Death array for UNK, TOT
## both TOT and UNK can be passed
deflateDeathsTOTUNKOAI <- function(Deaths){
  UNK <- TOT <- OPi <- OP  <- TOTTF <- UNKTF <- OPTF <- NULL
  
  TOTTF <- any(Deaths$Age == "TOT")
  if(TOTTF){
    TOT <- Deaths[Deaths$Age == "TOT", ]
    Deaths <- Deaths[Deaths$Age != "TOT", ]
  }
  UNKTF <- any(Deaths$Age == "UNK")
  if (UNKTF){
    UNK           <- Deaths[Deaths$Age == "UNK", ]
    Deaths        <- Deaths[Deaths$Age != "UNK", ]
  }  
  
  # slice off OP, to be slapped back on later
  OPi <- Deaths$AgeInterval == "+" 
  OPTF <- any(OPi)
  if (OPTF){
    OP            <- Deaths[OPi, ]
    Deaths        <- Deaths[!OPi, ]
  }
  return(
    list(Deaths=Deaths, UNK=UNK, OP=OP, TOT=TOT, UNKTF=UNKTF, OPTF=OPTF, TOTTF=TOTTF )
  )
}


library(HMDHFDplus)
CTYlist <- getHMDcountries()

for(iCTY in CTYlist){
  t.Deaths <- readInputDB(file.path("/data/wilmoth0/HMD/HMDWORK",iCTY), strict=FALSE)$Deaths
  t.Deaths <- deflateDeathsTOTUNKOAI(t.Deaths)$Deaths
  t.VV0 <- t.Deaths[t.Deaths$Lexis=="VV" & t.Deaths$AgeIntervali == 0, ]
  t.VV2p <- t.Deaths[t.Deaths$Lexis=="VV" & t.Deaths$AgeIntervali >1, ]
  cat("------------------")
  cat(iCTY,"\n")
  head(t.VV0)
  head(t.VV2p)
  cat("------------------")
}