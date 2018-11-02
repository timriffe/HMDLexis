
#' @title d_unk5 distributes differences between deaths1 and deaths2 to deaths1
#' all deaths are in standard internal format of the InDB

#' 
#' @param deaths1 death data frame  which will receive differenced death
#' @param deaths2 death data frame which contributes deaths to distribute proporitionally to deaths1
#' 
#' @return deaths data frame containing deaths1 adjusted for redistributed deaths from deaths2
#' 
#' @importFrom compiler cmpfun
#' @importFrom data.table 
#' 
#' @export
#' 

d_unk5 <- function(deaths1, deaths2){
  # verbatim translation of matlab script
  
  # convert to data.table with indexing by Year, Sex, AgeInterval
  
  minY=min(deaths1$Year)
  maxY=max(deaths1$Year)
  
  
  ## matlab original allocated wide age groups in d2 to d1, skipped UNK (-1), but allocated TOT (300) to TOT
  ## which is pointless.  The original skipped UNK deaths in d2 which is a mistake if the UNK category differs between 
  ## d2 and d2, so that category is not skipped.  The open age interval could differ between d1 and d2
  d.ans <- NULL
  for (i in  seq(from=minY, to=maxY)){ 
    cat(paste(i,"\n"))
    for(iSex in c("f","m") ){
      d1.in <- deaths1[ deaths1$Year==i & deaths1$Sex == iSex, ]
      d2.in <- deaths2[ deaths2$Year==i & deaths2$Sex == iSex, ]
      d1.out <- f.distribute_d2tod1(d1.in, d2.in)
      d.ans <- rbind(d.ans, d1.out)
    }
  }
  return(d.ans)
}
  
f.distribute_d2tod1 <- function(d1,d2){
  if(nrow(d2) == 0) return(d1)  # nothing to adjust
    
  for(j in seq(nrow(d2)) ){  # every entry in d2
      
    if(d2$Age[j] =="TOT"){
      next
    }
    if(d2$Age[j]=="UNK"){
      stopifnot( "UNK" %in% d1$Age) # nowhere to put UNK deaths
      if( d1[d1$Age=="UNK"] != d2$Deaths[j]) {
        message("d_unk5: unexpected UNK; ", d2$Deaths[j], " overwrites ", d1$Deaths[ d1$Age == "UNK"] )
      }
      d1[d1$Age=="UNK"] <- d2$Deaths[j]
      next
    }
    
    ## the below won't work if d1 age intervals overlap those of d2, e.g., [7,11) in d1 but [5,10) in d2
    ## in that case, an approach using cumulative death curves with interpolation might work
    if( d2$AgeInterval[j]=="+" || ( !is.na(d2$AgeIntervali)) ){  # open age interval or regular age entry
      isel1 <- !is.na(d1$Agei) & d1$Agei >= d2$Agei[j] & d1$Agei < d2$Agei[j] + ifelse(d2$AgeInterval[j]=="+" , 999, d2$AgeIntervali[j])
      if(any(isel1) ){          
        if( sum( d1$Deaths[isel1], na.rm=TRUE) != d2$Deaths[j]) {  # action needed to adjust death counts
          message("d_unk5: adjusting ", paste(d1$Deaths[isel1]), "deaths at ages ", paste(d1$Agei[isel1]), " by factor ", d2$Deaths[j]/ sum(d1$Deaths[isel1]))
          d1$Deaths[isel1] <- d1$Deaths[isel1] * (d2$Deaths[j]/ sum(d1$Deaths[isel1]))
        }
      }
      
    } else {
      warning("d_unk5: Expected a section of d1 in which to put d2 deaths, but not found")
    }
    
  } # end j
  
  return(d1)
  
} # f.distrubute_d2tod1

  # 
  # f.backcode<- function(x){  # old matlab coding for 
  #   x$MatAgei  <- x$Agei
  #   x$MatAgei <- ifelse( x$Age=="UNK", -1,  x$MatAgei)
  #   x$MatAgei <- ifelse( X$Age == "TOT", 300,x$MatAgei )
  #   
  # }