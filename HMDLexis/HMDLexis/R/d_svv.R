

#'@title a function to split deaths in vertical parallelograms.
#'
#'@description This method is described in an MP footnote, where we state that the method is to just split VV or VH in half. That's easy to implement, but the matlab code did much more. Since it was never vetted publicly, we do not implement the matlab oddities at this time, but perhaps a future maintainer will fill in that code chunk. Right now the code essentially essentially ignores the \code{reproduce.matlab}, argument, following the MP in all cases and returning a Warning message if \code{reproduce.matlab = TRUE}.
#' 
#' @param Deaths The standard Deaths data.frame, as typically returned by \code{readInputDB}
#'
#' @return Deaths with VV converted to TL and TU. There will be more rows now...
#' 
#' @export

d_svv <- function(Deaths, reproduce.matlab = FALSE){
  # this can be before or after redistributing deaths of Unknown age, it makes no difference
  # it can actually be at any time in the deaths processing (but before d_long(), obviosuly)
  
  # while we're at it, slice of UNK, rbind back on later:
  # using custom logical operator now with desired behavior
  DataKeep <- Deaths[Deaths$AgeIntervali %!=% 1 | Deaths$Lexis %!=% "VV", ]
  VV       <- Deaths[Deaths$AgeIntervali %==% 1 & Deaths$Lexis %==% "VV", ]
  if (reproduce.matlab){
    cat("\nWarning: you specified reproduce.matlab = TRUE.\nThe original matlab code did some out-of-protocol stuff that I didn't\nhave time to reproduce, sorry!\nThe function will return results for reproduce.matlab = FALSE.")
    reproduce.matlab <- FALSE
    
    for (sex in c("f","m")){
      
      # ignore infant TL
      TL <- DataKeep[DataKeep$Lexis %==% "TL" & DataKeep$Sex == sex & DataKeep$Agei %>% 0,] 
      TL$Agei <- TL$Agei + 1 # align to VV Age (it's confusing, but VV is indexed to the lower age)
      TU <- DataKeep[DataKeep$Lexis %==% "TU" & DataKeep$Sex == sex,]
      # get in AP matrices
      TL <- acast(TL, Agei ~ Year, sum, value.var = "Deaths", fill = NA_real_)
      TU <- acast(TU, Agei ~ Year, sum, value.var = "Deaths", fill = NA_real_)
      
      YrKeep  <- intersect(colnames(TU),colnames(TL))
      AgeKeep  <- intersect(rownames(TU),rownames(TL))
      # note that TU age remains, and TL age is -1 in case of VV alignment

    
      reproduce.matlab <- FALSE # provisional behavior until above is complete.
                                # looking for an efficient way to split using known proportions.
  }
  if (!reproduce.matlab){
    # in the MP, we say we divide equally into triangles..
    TL                <- TU              <- VV
    TL$Deaths         <- TL$Deaths / 2
    TU$Deaths         <- TU$Deaths / 2
    TL$Age            <- TL$Agei         <- TL$Agei + 1
    TL$Lexis          <- "TL"
    TU$Lexis          <- "TU"
    
    TLTU              <- rbind(TL,TU)
    TLTU$NoteCode1  	<- "d_svv()" # may overwrite previous NoteCode
  }
  
  D.Out <- d_agg(rbind(DataKeep, TLTU))
  
  invisible(resortDeaths(D.Out))
}















