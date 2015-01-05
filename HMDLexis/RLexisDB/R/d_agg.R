#' @title a function to aggregate late-registered deaths into the appropriate age and year categories
#' 
#' @description InputDB files typically contain a \code{"YearReg"} column, but it does not appear to be used either implicitly or explicitly in the LexisDB files. Further, oftentimes when a particular combination of Year-Sex-Age-Lexis is extracted from the data in the matlab code it is assumed that this will be a single number. In the case of additional late registered events, this operation will in fact return a vector, and could lead to unexpected results. This function is written to explicitly sum-in late registered events prior to further matlab operations, and can be considered a pre-processing step until such time as the matlab is phased out. In order to combine matlab and R in this way, we've decided to name the deaths file in the original raw state as \code{XYZdeath_pre.txt}, and the output of this function as the standard \code{XYZdeath.txt} so that the matlab programs will be able to find it. This is a temporary kludge for JPN and FIN, hopefully, although there are other countries where this procedure may be needed in the future, such as Finland or Australia, depending on how data come. This function may have other unexpected side-effects where redundant deaths have not been properly dealt with in the matlab code. For now, this function will be called on the CS-side, prior to submitting the InputDB, as this function needs to be used prior to completion of the LexisDB R programs.
#' 
#' @param Deaths.In This should be the standard Deaths data.frame, as typically returned by \code{readInputDB()}. 
#' 
#' @importFrom data.table data.table
#' 
#' @export 
#' 
#Deaths.In <- read.table("/data/commons/triffe/git/HMD_CS/HMDwork/C_JPN/JPN/InputDB/JPNdeath.txt", sep = ",", na.strings = ".", header = TRUE, stringsAsFactors = FALSE)
## take only LDB = 1
#Deaths.In <- Deaths.In[Deaths.In$LDB == 1, ]
## don't deal with UNK year yet (that's easy, will do next)
#Deaths.In <- Deaths.In[Deaths.In$Year != "UNK", ]
#Deaths.In <- Deaths.In[!is.na(Deaths.In$YearInterval), ]
#
## don't deal with open year yet (less easy, will put this off until further discussion)
#Deaths.In <- Deaths.In[Deaths.In$YearInterval != "-", ]
d_agg <- function(Deaths.In){
  # remove TOT (never needed)
  Deaths.In <- Deaths.In[Deaths.In$Age != "TOT", ]
  
  # aggregate by all potentially relevant columns
  D.Agg            <-  data.table(Deaths.In)[, sum(Deaths), 
                          by = c("Year","YearInterval","Sex","Age","AgeInterval","Lexis")]
            # library(data.table)  
  
  setnames(D.Agg,"V1","Deaths")
  head(D.Agg[D.Agg$Year == 2001 & D.Agg$Age == "0", ],11)
  head(Deaths.In[Deaths.In$Year == 2001 & Deaths.In$Age == "0" & Deaths.In$YearReg == 2003, ],15)
  
#  D.Agg            <-  data.table(Deaths.In)[, function(.SD){
#      if (nrow(.SD) <= 1){
#        return(eval(copy(.SD)))
#      } else {
#        out <- copy(.SD[1, ])
#        out[["Deaths"]]  <- sum(.SD[["Deaths"]])
#        out[["RefCode"]] <- paste(.SD[["RefCode"]], sep = ",")
#        if (any(!is.na(.SD[["NoteCode1"]]))){
#          out[["NoteCode1"]] <- paste(.SD[["NoteCode1"]][!is.na(.SD[["NoteCode1"]])],sep=",")
#        }
#        if (any(!is.na(.SD$NoteCode2))){
#          out[["NoteCode2"]] <- paste(.SD[["NoteCode2"]][!is.na(.SD[["NoteCode2"]])],sep=",")
#        }
#        out[["NoteCode3"]] <- "d_agg()"
#        out
#    }
#  }, 
#                          by = c("Year","YearInterval","Sex","Age","AgeInterval","Lexis")]
#   D.Agg            <-  data.table(Deaths.In)[, .I[Deaths=sum(Deaths)], 
#                          by = c("Year","YearInterval","Sex","Age","AgeInterval","Lexis")]
#                        
#   test <- data.table(a=c(1,2,3,4,5),b=letters[1:5],g=c(1,1,3,3,3))
#   test[,.I[sum(a)],by=g]                
#            ?data.table            
#  D.agg <-    do.call(rbind,
# lapply(
#   split(Deaths.In, f=with(Deaths.In, list(Year,YearInterval,Sex,Age,AgeInterval,Lexis))),
#   function(DAT){
#     if (nrow(DAT) == 1){
#       return(DAT)
#     } else {
#       out <- DAT[1,]
#       out[["Deaths"]]  <- sum(DAT[["Deaths"]])
#       out[["RefCode"]] <- paste(DAT[["RefCode"]], sep = ",")
#       if (any(!is.na(DAT[["NoteCode1"]]))){
#         out[["NoteCode1"]] <- paste(DAT[["NoteCode1"]][!is.na(DAT[["NoteCode1"]])],sep=",")
#       }
#       if (any(!is.na(DAT$NoteCode2))){
#         out[["NoteCode2"]] <- paste(DAT[["NoteCode2"]][!is.na(DAT[["NoteCode2"]])],sep=",")
#       }
#       out[["NoteCode3"]] <- "d_agg() has been run on these data!"
#       out
#     }
#   }))
                        
  # V1 is Deaths
  #colnames(D.Agg)[colnames(D.Agg) == "V1"] <- "Deaths"

  # fill in unused columns with NA
#  D.Agg$NoteCode3  <- D.Agg$NoteCode2 <- D.Agg$RefCode <- 
#    D.Agg$Access   <- D.Agg$YearReg <- NA
#  
#  # make note that this is no longer an input file
#  D.Agg$NoteCode1  <- "d_agg() has been run on these data!"

  # resort columns
  D.Out            <- D.Agg[,colnames(Deaths.In)]
  
  # return invisibly
  invisible(D.Out)
}
