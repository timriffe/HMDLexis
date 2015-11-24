
library(magrittr)
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

# process populations:

tmp <- Pop
Pop <- tmp
Pop     <- p_unk(Pop)
Pop<-   p_ic(Pop, Deaths, Births, reproduce.matlab = TRUE)

Pop<-   p_srecm(Pop, Deaths, reproduce.matlab = TRUE)



# generate LDB output:

# need writeLDB() function here.
LDB.struc <- computeLDB(Births, Deaths, Pop, 
                        Tadj= (if(TadjTF){ Tadj} else NULL), verbose=FALSE) 

writeLDB(LDB.struc,  LDBpath=file.path(xWORKING,"LexisDB"), COUNTRY=LDB.struc$COUNTRY, strict=TRUE)








