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
#' @importFrom data.table rbindlist
#     note that rbindlist() is the same as do.call(rbind, my.list), which comes in base
#     rbindlist() does not check column names, FYI. Maybe we don't care about that.

# TR: I've made no changes to the function Feb 2, 2016

## WORKING = "/home/boe/DEMOG/Triffe_git/HMDLexis/HMDLexis"

require(data.table)
computeLDB <- function(  # WORKING = "/data/commons/hmd/HMDWORK/DNK",
                        Births,
                        Deaths,
                        Pop,
                        Tadj,     # Tadj is required, even if all defaults of 1 are used
                        #save.bin = FALSE,
                        verbose = FALSE){

  ## processes a single country.  Test that parts agree
  Tadj.name <- as.character(Tadj$PopName)
  
  if( length(unique( c( Pop$PopName, Deaths$PopName, Births$PopName, Tadj.name))) != 1) {
    stop("More than one PopName detected among components")
  }
  
  COUNTRY <- unique(Pop$PopName)

  AGES.min <- 0
       ## when 'strict', max age is 130, but can stretch to 150 in some pops like USA
       ## exclude Tadj from this computation, since the ages there are a result of computations
       ## and do not necessarily conform
  AGES.max <- max(Pop$Agei, Deaths$Agei, na.rm=TRUE )
  AGES <- (AGES.min):(AGES.max)
  MaximumAge <- AGES.max
  ## Tadj$Age needs to be a superset of AGES
  print("length( setdiff(AGES, as.numeric(Tadj$Age )) ) == 0") 
  print( length( setdiff(AGES, as.numeric(Tadj$Age )) ) )
  stopifnot( length( setdiff(AGES, as.numeric(Tadj$Age )) ) == 0)
  
  SEXES <- c("f","m")

  FirstYear <- min(Pop$Year, Deaths$Year)
  LastYear  <- max (Pop$Year,Deaths$Year)
  YEARS  <- seq(from=FirstYear, to=LastYear, by=1)
  ## create 0:130 template of output, then cycle through years, sex and slot the values into place
  ## special process for births.  Append each year x sex to output array
  ## Year, Age, Triangle, Cohort, Populaton, Deaths
  this.template<- data.frame(Year=2000+NA, Age=AGES,Triangle=1+NA, Cohort=2000+NA,
                             Population=0+NA, Deaths=0+NA )

 

  result<- this.template[0,]
  
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
        this.Pop <- this.Pop[ order(this.Pop$Agei),]

        if( any(this.Pop$Age %in% c("UNK","TOT")) ){
          stop(paste("UNK or TOT left in Pop array, year=", year))
        }

        stopifnot( testAgeMatch(this.Pop$Agei, AGES, "Pop", year, sex) )

        ##  will have zero rows if year out of range
        this.Pop.next <- Pop[ Pop$Sex==sex & Pop$Year == year+1 & Pop$LDB==1, ]
        if(nrow(this.Pop.next)==0){  ## create fake Pnext if out of data
          this.Pop.next <- this.Pop
          this.Pop.next$Year <- this.Pop.next$Year + 1
          this.Pop.next$Population <- rep( -1L, nrow(this.Pop.next) )
        }
        stopifnot( testAgeMatch(this.Pop.next$Agei, AGES, "Pop.next", year, sex) )
        
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
        
        this.DeathsU.full.minus.1 <-  rbind(c(-0L, NA),this.DeathsU.full[-1,] )
        this.DeathsU.full.minus.1$Agei <- this.DeathsU.full.minus.1$Agei -1L
        
       
        this.B <- Births[ Births$Year == year & Births$Sex == sex & Births$LDB == 1, c("Area","Births")]
        if( NROW(this.B) >= 2){
          stop(paste("births multiply defined, year, sex:", year, sex))
        }
        if( NROW(this.B) != 1 ){
          warning(paste("births missing or multiply defined for year, sex", year, sex, this.B, "\n using -1"))
          this.B <- data.frame(Area="", Births = -1)
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
          this.Vx <- Tadj[ Tadj$Type=="Vx" & Tadj$Year == year & Tadj$Sex == sex & 
            Tadj$LDB == 1, c("Age", "Area1", "Value")]
          Vx <-  rep(1L, length(this.Pop$Agei))
          Vx[ match(this.Vx$Age, this.Pop$Agei, nomatch=0)] <- this.Vx$Value[ match(this.Vx$Age ,  this.Pop$Agei, nomatch=0) ]
        
          
          ##  TODO: test for consistency / multiplicity of Area values, but ignore for now
          
        }
        ## exposures (lives crossing indexed age during year; age 0 handled separately)
        ## N=(p(indp,4)-du+pnext(indpn,4)/Vx+dl)/2
        
        this.N <- ( ( this.Pop$Population + this.DeathsL.full$Deaths) +
                    ( this.Pop.next$Population / this.Vx$Value  - this.DeathsU.full.minus.1$Deaths)  ) / 2
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
    
    result.by.sex[[isex]] <- data.table::rbindlist(result.by.year)  # much faster version of rbind() for data.frames
    
    
  } # end isex

  ## return structure
  return(list( f=result.by.sex[["f"]], m=result.by.sex[["m"]], COUNTRY=COUNTRY) )
} #end computeLDB





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
