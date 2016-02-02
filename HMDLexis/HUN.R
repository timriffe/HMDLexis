
#library(magrittr)
getwd()
devtools::load_all("HMDLexis")

PROJECT = "HMDWORK"  # top level project directory
XXX <- "HUN"
xWORKING = file.path(PROJECT,XXX)
## for development
source("computeLDB.R")
source("writeLDB.R")


# read in InputDB
IDB     <- readInputDB(WORKING = xWORKING,
                  XXX = NULL,
                  log.file = "", # optional logging of input errors
                  InputDB.name = "InputDB",
                  save.bin = FALSE, verbose = TRUE)

# get individual objects
Births  <- IDB$Births
Deaths  <- IDB$Deaths
Pop     <- IDB$Pop

# TR: readInputDB() always produces a full tadj, containing 1s for Vx
# for all years that have no change, and using the input Vx values where
# required. So this conditional not required. It's safe to assign a Tadj
TadjTF  <- FALSE
if(! is.null(IDB$Tadj)){
  Tadj <- IDB$Tadj
  TadjTF <- TRUE
}

## > head(Tadj,3)
##   PopName Year Age Area1 Area2 Sex Type Value LDB
## 1     DNK 1835  NA     2     2   f   Rb     1   1
## 2     DNK 1835   0     2     2   f   Vx     1   1
## 3     DNK 1835   1     2     2   f   Vx     1   1
# This is the nested, less legible way:
# Deaths <- d_long(d_soainew(d_unk(Deaths)))
# same as below in magrittr (using piping):

tmp<- Deaths
Deaths <- d_unk(Deaths)
Deaths <- d_soainew(Deaths)
Deaths <- d_long(Deaths)
any(is.na(Deaths$Deaths))
#length(unique(Deaths$Year)) * 
#  length(unique(Deaths$Sex)) * 
#  length(unique(Deaths$Age)) * 
#  length(unique(Deaths$Lexis)) == nrow(Deaths)
#length(unique(Deaths$Year)) * 2 * 131 * 2 == nrow(Deaths)
# process populations:

tmp <- Pop
Pop <- tmp
Pop     <- p_unk(Pop)
Pop<-   p_ic(Pop, Deaths, Births, reproduce.matlab = TRUE)

Pop<-   p_srecm(Pop, Deaths, reproduce.matlab = TRUE)

#length(unique(Pop$Year)) * 
#  length(unique(Pop$Sex)) * 
# length(unique(Pop$Age))  == nrow(Pop)
#length(unique(Pop$Year)) * 2 * 131 == nrow(Pop)
#unique(Pop$Area)
#head(Pop[is.na(Pop$Area),])
# generate LDB output:
# TR: could do Pop <- p_long(Pop) redundantly in this case
# need writeLDB() function here.
# TR: see comments in header of computLDB()
LDB.struc <- computeLDB(Births, Deaths, Pop, 
                        Tadj= (if(TadjTF){ Tadj} else NULL), verbose=FALSE) 

writeLDB(LDB.struc,  LDBpath=file.path(xWORKING,"LexisDB"), COUNTRY=LDB.struc$COUNTRY, strict=TRUE)








