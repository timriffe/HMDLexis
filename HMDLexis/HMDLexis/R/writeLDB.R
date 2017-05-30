## TODO:  allow option to supress header of CSV file (the old format)
###############################################################################
writeLDB <- function(ldb.by.sex, 
                     LDBpath,  # WORKING/XXX/LexisDB   typically
                     COUNTRY=NULL,
                     #save.bin = FALSE,
                     strict = TRUE){
 
  dir.create(LDBpath,recursive=TRUE, mode="0775")
  
  # MatlabRound() is for rounding output, should give same result as matlab, assuming that's important
  # by CB, updated by TR to take digits as arg.
  # CAB.  Fixed rounding formula to handle negative numbers
  MatlabRoundLDB <- function(x, digits = 0){ 
    
    # this 1) rounds, and
    # 2) makes sure that the zeros stay on the end of the number string to the specified number of digits
    
    ## TODO: use tapply to ensure that any rounding sums to total?
    if (is.numeric(x)){
      fac <- rep(10 ^ digits, length(x))
      #x   <-  floor(x * fac + sign(x) * 0.5) / fac
      x   <-  sign(x) * ( floor( abs(x) * fac + 0.5)/ fac )
    }
    return(x)
  }
  
  if(is.null(COUNTRY))
    COUNTRY<- ldb.by.sex$COUNTRY
  
  ldb.f <- ldb.by.sex[["f"]]
  ldb.m <- ldb.by.sex[["m"]]
  
  if(strict){
    ldb.f$Population <- MatlabRoundLDB(ldb.f$Population, digits=2)
    ldb.f$Deaths <- MatlabRoundLDB(ldb.f$Deaths, digits=2)
    ldb.m$Population <- MatlabRoundLDB(ldb.m$Population, digits=2)
    ldb.m$Deaths <- MatlabRoundLDB(ldb.m$Deaths, digits=2)
  } 
    
  write.csv(ldb.f, 
              file=file.path(LDBpath, paste0("f",COUNTRY,".txt") ),
              eol="\r\n",
              row.names=FALSE)
  write.csv(ldb.m, 
              file=file.path(LDBpath, paste0("m",COUNTRY,".txt") ),
              eol="\r\n",
              row.names=FALSE)    
    
  }
  
