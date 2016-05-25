


#' @title compute the standard LDB object for males and females
#' @description the LDB object has columns for Year, Age, Lexis (1,2), Population, and Deaths. Lexis = 2 means upper triangle deaths and Jan 1 Populations. Lexis = 1 means lower triangle deaths and lower age bound birthdays (the bottom horizontal line of the Lexis square). This is one way to capture births for Age = 0, but Lexis = 1 Population counts are otherwise not used in any calculations. TR has proposed to instead put DL, DU in columns, and Pop Jan 1, Pop Dec 31 in columns, which would obviate the need for Tadj from the InputDB downstream. Further, new exposure formulas count be accomodated by sd and bbar columns, rather than calculating on the fly, then no need to read in Monthly births from InputDB for LifeTables calcs. In either case, all R \code{data.frame} arguments should be fully processed before getting to this stage.
#' 
#' @param Births the standard R Births \code{data.frame}. All anomalies should be accounted for before this function.
#' @param Deaths the standard R Deaths \code{data.frame}, fully processed and complete.
#' @param Pop the standard R Population \code{data.frame}, fully processed and complete.
#' @param Tadj the standard R Territorial Adjustment \code{data.frame}, as returned by \code{readInputDB()}.
#' @param verbose. logical. do you want lots of info to go to the console?
#' 
#' @export
#' @importFrom data.table rbindlist
#' 

computeLDB <- function(  # WORKING = "/data/commons/hmd/HMDWORK/DNK",
  Births,
  Deaths,
  Pop,
  Tadj,     
  verbose = FALSE){
  
  ## processes a single country.  Test that parts agree
  if(TadjTF){
    Tadj.name <- as.character(Tadj$PopName)
  } else {
    Tadj.name <- NULL
  }
  
  if( length(unique( c( Pop$PopName, Deaths$PopName, Births$PopName, Tadj.name))) != 1) {
    stop("More than one PopName detected among components")
  }
  
  COUNTRY     <- unique(Pop$PopName)
  
  AGES.min    <- 0
  ## when 'strict', max age is 130, but can stretch to 150 in some pops like USA
  AGES.max    <- max(Pop$Agei, Deaths$Agei, as.numeric(Tadj$Age ), na.rm=TRUE )
  AGES        <- (AGES.min):( min(AGES.max, 130) )
  MaximumAge  <- AGES.max
  SEXES       <- c("f","m")
  
  FirstYear   <- min(Pop$Year, Deaths$Year)
  LastYear    <- max (Pop$Year,Deaths$Year)
  YEARS       <- seq(from=FirstYear, to=LastYear, by=1)
  ## create 0:130 template of output, then cycle through years, sex and slot the values into place
  ## special process for births.  Append each year x sex to output array
  ## Year, Age, Triangle, Cohort, Populaton, Deaths
  this.template<- data.frame(Year=2000+NA, Age=AGES,Triangle=1+NA, Cohort=2000+NA,
    Population=0+NA, Deaths=0+NA )
  
  result <- this.template[0, ]
  
  result.by.sex <- vector("list", length=length(SEXES))
  names(result.by.sex) <- SEXES
  
  for( isex in seq(along=SEXES)){
    sex <- SEXES[isex]
    
    result.by.year <- vector("list", length=length(YEARS))  # hold results
    names(result.by.year) <- YEARS
    
    for( iyear in seq(along=YEARS) ){
      
      year<- YEARS[iyear]
      thisLower <- this.template
      thisUpper <- this.template
      
      this.Pop <- Pop[ Pop$Sex==sex & Pop$Year == year & Pop$LDB==1, ]
      ## CDW - population vector needs to be explicitly ordered by age
      this.Pop <- this.Pop[ order(this.Pop$Agei),]
      
      if( any(this.Pop$Age %in% c("UNK","TOT")) ){
        stop(paste("UNK or TOT left in Pop array, year=", year))
      }
      
      stopifnot( testAgeMatch(this.Pop$Agei, AGES, "Pop", year, sex) )
      
      ##  will have zero rows if year out of range
      this.Pop.next <- Pop[ Pop$Sex==sex & Pop$Year == year+1 & Pop$LDB==1, ] # CB, but rows do not have ages equal to ages of other arrays, shifted by one.  Awkward.
      ## CW adds this next line to explicitly sort this.Pop.next by age 
      this.Pop.next <- this.Pop.next[ order(this.Pop.next$Agei),]
      if(nrow(this.Pop.next)==0){  ## create fake Pnext if out of data
        this.Pop.next <- this.Pop
        this.Pop.next$Year <- this.Pop.next$Year + 1
        this.Pop.next$Population <- rep( -1L, nrow(this.Pop.next) )
      }
      stopifnot( testAgeMatch(this.Pop.next$Agei, AGES, "Pop.next", year, sex) )
      ## CDW - compute population vector, offset by one year of age
      ## Agei[1] = -1 (dummy), Agei[2] = 0, Agei[3] = 1, etc....
      this.Pop.minus.1 <-  rbind(c(-0L, NA),this.Pop[-1,] )
      this.Pop.minus.1$Agei <- this.Pop.minus.1$Agei -1L
      
      #d=selif(deaths, deaths(:,1)==i & deaths(:,3)==sex);
      this.DeathsL <- Deaths[ Deaths$Year==year & Deaths$Sex == sex & 
          Deaths$Lexis=="TL" & Deaths$LDB==1, ]
      this.DeathsU <- Deaths[ Deaths$Year==year & Deaths$Sex == sex & 
          Deaths$Lexis=="TU" & Deaths$LDB==1, ]
      
      ## Deaths can have NA above a certain age...not necessarily conformable 
      ##stopifnot( testAgeMatch(this.DeathsL$Agei, AGES, "DeathsL", year, sex) )
      ## TODO: test for two entries for each age (TL, TU)
      
      if(verbose){
        testAgeMatch(this.DeathsL$Agei, AGES, "DeathsL", year, sex) 
        testAgeMatch(this.DeathsU$Agei, AGES, "DeathsU", year, sex) 
      }
      
      ## this will replace in Y the rows where X has matching t
      # Y[  match(X$t, Y$t, nomatch=0), ] <-  X [ X$t %in% Y$t,]
      
      ## warning: a bit tricky.  Done this way to avoid possible problems with conformability
      ## and sorting of ages
      this.DeathsL.full <- data.frame(Agei=AGES, Deaths=NA)  # need full age range, conformable to Pop
      this.DeathsL.full[ match(this.DeathsL$Agei, this.DeathsL.full$Agei, nomatch=0), ] <-
        this.DeathsL[ this.DeathsL$Agei %in% this.DeathsL.full$Agei,c("Agei", "Deaths")]
      
      this.DeathsU.full <- data.frame(Agei=AGES, Deaths=NA)  # need full age range, conformable to Pop
      this.DeathsU.full[ match(this.DeathsU$Agei, this.DeathsU.full$Agei, nomatch=0), ] <-
        this.DeathsU[ this.DeathsU$Agei %in% this.DeathsU.full$Agei,c("Agei", "Deaths")]
      
      # CB: this calculaton is wrong! Unless we match on ages in some way.
      # are we assuming equivalent ages are in the same row across objects?
      this.DeathsU.full.minus.1 <-  rbind(c(-0L, NA),this.DeathsU.full[-1,] )
      this.DeathsU.full.minus.1$Agei <- this.DeathsU.full.minus.1$Agei -1L
      
      
      this.B <- Births[ Births$Year == year & Births$Sex == sex & Births$LDB == 1, c("Area","Births")]
      if( NROW(this.B) >= 2){
        stop(paste("births multiply defined, year, sex:", year, sex))
      }
      if( NROW(this.B) != 1 ){
        warning(paste("births missing or multiply defined for year, sex", year, sex, this.B, "\n using -1"))
        this.B <- data.frame(Area=NA, Births = -1)
      }
      
      ## if Area designations vary within year, all bets are off
      this.Areas <- unique(c(this.B$Area, this.Pop$Area, this.DeathsU$Area, this.DeathsL$Area))
      this.Areas <- this.Areas[ ! is.na(this.Areas)]
      if( length(this.Areas) != 1){
        ## p_ecm() does not yet preserve areas, so warn rather than abort
        warning("multiple Area designations in year, sex Pop, Deaths, Births data")
      }
      
      ## Tadj adjustment factors.  Use default value of 1 unless overridden
      ##  GenericTadj is handled when InputDB values are read in, so this does not 
      ## need repeating
      
      if(is.null(Tadj)){
        ## defaults
        this.Rb <- 1
        this.Vx <- cbind(this.Pop$Agei,1)
        warning("empty Tadj encountered")
        
      } else {
        stopifnot( ! any(Tadj$Type=="Rd"))  # not really sure what this is or how to handle it.  Deaths seem to be always reported as adjusted.  Calcs involving deaths do not span multiple years, so never trigger use of Tadj factors
        
        this.Rb <- Tadj[Tadj$Type=="Rb" & Tadj$Year == year & Tadj$Sex == sex & 
            Tadj$LDB == 1, c("Area1", "Value")]
        
        ## CAB: test this---I am getting length 151 for this.Vx; why do Tadj ages span 0:150 ?
        # TR (25 May, 2016): because the USA had a Tadj file that for some crazy reason went to age 150
        # readInputDB() by default would throw an error if there is anything greater than 130.
        # however, since the USA was a test case and we wanted to actually be able to test
        # I put in an argument called 'strict' (logical). if TRUE, then we enforce no ages greater than 130
        # If FALSE we let Tadj go up to 150 to account for this silly legacy without breaking
        this.Vx <- Tadj[ Tadj$Type=="Vx" & Tadj$Year == year & Tadj$Sex == sex & 
            Tadj$LDB == 1, c("Age", "Area1", "Value")]
        Vx <-  rep(1L, length(this.Pop$Agei))
        Vx[ match(this.Vx$Age, this.Pop$Agei, nomatch=0)] <- this.Vx$Value[ this.Vx$Age %in% this.Pop$Agei ]
        #^^^  problem here with NAs from match
        
        ##  TODO: test for consistency / multiplicity of Area values, but ignore for now
        
      }
#------------------------------------------------------------      
# TR (25 May, 2016): are these really exposures? I thought it was calculating the 
# birthday line, a la:
#
#      |--------|
#      |     /  |
#      |   /    |P2
#      | /   DL |
#      |---Nx---| <- Nx birthday line
#      | DU  /  |
#    P1|   /    |
#      | /      |
#      |--------|
# Nx = ((P1 - DU) + (P2 + DL)) / 2   
# from the looks of the code below something is swapped,
# but I'm not sure really, and I never did look into 
# Lexis = 1 Population counts, since they aren't used per se.
# although there exists a redux of the v5 exposure formula
# that would give an identical result.
#------------------------------------------------------------

      ## exposures (lives crossing indexed age during year; age 0 handled separately)
      ## N=(p(indp,4)-du+pnext(indpn,4)/Vx+dl)/2
      ## CDW - update formula above: p(indp,4) should be this.Pop.minus.1
      ## not this.Pop
      this.N <- ( ( this.Pop.minus.1$Population + this.DeathsL.full$Deaths) +
          ( this.Pop.next$Population / Vx  - this.DeathsU.full.minus.1$Deaths)  ) / 2
      ##indb=find(b(:,2)==sex & b(:,3)==i);
      ##if isempty(indb)
      ##N=-1;
      ##else
      ##  N=b(indb,4);
      ## end
      
      
      this.N[1] <- this.B$Births   # age 0
      
      # Year, Agei, Cohort=Year-Agei, N, Dl
      #thisLower$Year <-  year
      #thisLower$Age  <-  AGES
      #thisLower$Cohort <- thisLower$Year - thisLower$Age
      
      ##    indp=find(p(:,2)==0);
      ##    if isempty(indp)
      ##    P=-1;
      ##    else
      ##      P=p(indp,4);
      ##    end
      
      ##fprintf(fout, '%g, %g, 1, %g, %.2f, %.2f\n', i, j, i-j, N, Dl);
      ##fprintf(fout, '%g, %g, 2, %g, %.2f, %.2f\n', i, j, i-j-1, P, Du);
      
      this.DT.1 <- data.frame( Year=year, Age=this.Pop$Agei, Triangle = 1, 
        Cohort = year - this.Pop$Agei, Population = this.N, Deaths = this.DeathsL.full$Deaths)
      this.DT.2 <- data.frame( Year=year, Age=this.Pop$Agei, Triangle = 2, 
        Cohort = (year - this.Pop$Agei - 1), Population = this.Pop$Population, 
        Deaths = this.DeathsU.full$Deaths)
      this.DT <- rbind(this.DT.1 , this.DT.2)
      this.DT <- this.DT[ order(this.DT$Age, this.DT$Triangle), ]
      ## replace NAs with -1 for value cols
      this.DT$Population <- ifelse( is.na(this.DT$Population), -1, this.DT$Population)
      this.DT$Deaths <- ifelse( is.na(this.DT$Deaths), -1, this.DT$Deaths)
      
      
      ## store results
      result.by.year[[iyear]] <- this.DT
      
    }  # end iyear
    #browser()
    result.by.sex[[isex]] <- rbindlist(result.by.year)  # much faster version of rbind() for data.frames
    
    
  } # end isex
  
  ## return structure
  return(list( f=result.by.sex[["f"]], m=result.by.sex[["m"]], COUNTRY=COUNTRY) )
} #end computeLDB


#' @title an age checker
#' @description this is a utility used many times in \code{computeLDB()}
#' @param ages1 first set of ages (lower bounds)
#' @param ages2 second set of ages (lower bounds)
#' @param desc optional, set of the fly to identify errors
#' @param tyear optional, set of the fly to identify errors
#' @param tsex optional, set of the fly to identify errors
#' 
#' @export
testAgeMatch<- function(ages1, ages2, desc, tyear=NULL, tsex=NULL){
  # assert that ages match
  res<- TRUE
  if(length(ages1) != length(ages2)){
    cat(paste("Age classes lengths not equal", desc, "year, sex ", tyear, tsex, "\n"))
    res<- FALSE
  }
  if(length( setdiff(ages1, ages2) != 0 )){
    warning(paste( "age mismatch in data for ", desc, " year, sex", tyear, tsex) )
    cat(paste("Age classes mismatch:", desc, "year, sex ", tyear, tsex, "\n"))
    cat( paste("Ages mismatched: ", paste(setdiff(ages1, ages2), collapse=""), "\n" ) )
    res<- FALSE
  }
  return(res)
}

## Q's  
## (1) do Tadj  ever really enter?   Are they not used prior to LexisDB output?  Seems bad
##     to mix adjustment of data with output of data
# TR:       I think yes, because the Population column of pop only refers to January 1 by the end of
#           of processing. For this reason, you still need Tadj to calculate. I wish that LDB column for Lexis
#           (1 or 2) referred to Jan 1 (1), Dec 31 (2) instead of Jan 1 & AP birthday line... We only
#           gain an entry for births by this decision, but they are not used in Lifetable calcs. Jan 1 and Dec 31
#           are however used. If we made this change then the LDB file would be sufficient for v5 calcs.
#           otherwise the Lifetable functions need to dip into InputDB to tadj for Dec 31...

## (2) CB: confirm with TR, but it appears that there is no need to test for TOT and UNK cases
##     since those lines are purged before this ultimate output step
# TR:       most pop functions redundantly remove TOT. Only p_unk() distributes UNK, so checking for UNK
#           would only tell you whether it needs to be run. In that case, you could just call it, and in-trickles
#           the dream of a fully automatic LDB machine. Turing capable and everything.

## (3) Deaths have NA above age 99.  Like #1, surprising that data is not manipulated before
##     i/o outputting the data, and that the age ranges differ for Pop and Deaths   Is NA used
##     for missing in the case of Pop?
# TR:       HUN does not produce this. This must be a bug by some production function. Gotta find it.

## (4) type "Rd" is a possibility in Tadj, yet this is not readily accounted for in the codebase, and
##     only occurs for a single country.  
# TR:       correct, is it even in the MP? If not, add to list for v7.

## (5) Area consistency -- not really enforced.  Assumption that Tadj$Area1 will match up with Area fields
##     in Births, Pop, Deaths, but this is not enforced much.  Prior routines (p_ic() ?) augment the original data
##     but do not nicely track Area fields, with the result that the Pop data contains many NAs.  This starts at
##     age 80 and increases in frequency with higher ages.  Proper incorporation of Tadj remains TODO in many 
##     functions and it is too late by this stage to make good use of Tadj.

# TR:       see new function assignArea() in utils.R  This has been incorporated into p_ecm() and p_srm()
#           So HUN now correctly computes in this respect. Checking consistency between inputs needs to happen
#           in readInputDB(). I've not yet done this. Feel free.

## (6) p_ecm() routine does not preserve 'Area', so change error to warning for now
# TR: fixed

##
## (7) split function of matlab routines into 2 parts.  First is computation and second is I/O.
##     NB: Matlab rouding algorithm needed correction, as it did not handle -1 correctly (or any negative #s)
# TR: NB: one common difference that is toggled by reproduce.matlab is early rounding. Further, now that we have
#      changes in exposure and a0 calcs, nearly every Lifetable entrey will be slightly different, and it's the 
#      perfect time to depart from Matlab rounding norms, which we've not fully reproduced anyway. Why try? 

##writeLDB()

## uses routine from data.table package, assumed loaded
# TR: add roxygen header, like this:
##' @importFrom data.table rbindlist
#     note that rbindlist() is the same as do.call(rbind, my.list), which comes in base
#     rbindlist() does not check column names, FYI. Maybe we don't care about that.

# TR: I've made no changes to the function Feb 2, 2016

## WORKING = "/home/boe/DEMOG/Triffe_git/HMDLexis/HMDLexis"

## CW - April 21, 2016 - edits made to population inputs and to the computation of N:  

