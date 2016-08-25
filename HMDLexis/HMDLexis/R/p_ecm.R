#'
#' @title p_ecm a simple extinct cohort method function, for modularity
#' 
#' @description This implementation seeks to follow the MP. Deaths must be finalized. Populations must be in single ages and all years (i.e. after \code{p_split()} and \code{p_ic()}, if necessary). This function estimates population counts from death counts within area B from Figure 6 'Methods used for population estimates'. The actual work is odone by \code{p_ecm_inner()}. This function can handle territorial or universe adjustments. Results match those of matlab exactly, except at the time of implementation, for cohorts passing through a territorial adjustment in years prior to the adjustment things are slightly different due to a temporary bug in matlab, which will probably be fixed in due time. This version has been hand-checked to do the right thing when \code{Tadj} is specified.
#' 
#' @param Pop The standard internal Population data.frame, as
#' @param Deaths after all processing is done. Completed triangles.
#' @param Tadj LexisDB internal territorial (or universe) adjustment \code{data.frame}, as given by \code{readInputDB()}. This is optional.
#' @param a lowest age for which EC estimates should be given. Default 80.
#' @param reproduce.matlab logical. default is \code{FALSE}. This affects only the border cohort between EC and SR.
#' 
#' @return Population with old ages either filled in or imputed using the extinct cohort method
#' 
#' @details Conceivably, a population could have a long series with low open ages in the early part, but decent population data in the recent part. If we want a value of \code{A} lower than 80 for only part of the series, this can be achieved by subsetting \code{Deaths} for different year ranges, and making two calls to the function, then \code{rbind()}ing back together, removing duplicates as necessary. \code{reproduce.matlab} only affects the starting value of the border cohort between EC and SR methods, not a really big deal. Setting to \code{TRUE} does not activate the bug known from matlab.
#' 
#' @export
#' 
 
p_ecm <- function(Pop, Deaths, Tadj = NULL, a = 80, reproduce.matlab = FALSE){
	
	ColnamesKeep 	<- colnames(Pop)
	
	Pop 			<- Pop[Pop$Age != "TOT", ]
	
	# TR: new Aug 25, 2016 for Tadj
	if (! is.null(Tadj)){
		# only keep Vx for Tadj
		Tadj 	      <- Tadj[Tadj$Type == "Vx", ] 
		Tadj$Cohort   <- Tadj$Year - Tadj$Age - 1
	}
	
	# slice off UNK, rbind back on later:
	UNKi          <- Pop$Age == "UNK"
	UNKTF         <- any(UNKi)
	if (UNKTF){
		UNK           <- Pop[UNKi, ]
		Pop           <- Pop[!UNKi, ]
	}
	
	# add Cohort Columns for alignment
	# MP uses omega, but we use the most recent extinct cohort, since
	# we mark Cohorts and can use them to select.
	Deaths        <- d_addCohortColumn(Deaths)
	Pop           <- p_addCohortColumn(Pop)
	
	# Sex <- "f"; Sex <- "m"
	PopMF <- list()
	for (Sex in c("m","f")){
		Dsex              <- Deaths[Deaths$Sex == Sex, ]
		Psex              <- Pop[Pop$Sex == Sex, ]
		
		# TR: new Aug 25, 2016 for Tadj
		if (! is.null(Tadj)){
			Tsex              <- Tadj[Tadj$Sex == Sex, ]	
		} else {
			Tsex              <- Tadj
		}
		
		# 1) determine which cohorts are extinct. 
		#For now assume that this is separate for each sex
		omega             <- p_ecm_findOmega(Dsex, l = 5, threshold = 0.5)
		# p_ecm_inner() defined below, in same script. not quirky, but it is useful
		# to have these steps be modular, so that it can be called elsewhere
		ECpop             <- p_ecm_inner(
				Dsex = Dsex, 
				Tsex = Tsex, 
				a = a, 
				omega = omega, 
				reproduce.matlab = reproduce.matlab)
		
		# TR: p_ecm() was guilty of not assigning Area. Since future Tadj compatibility
		# will mostly be accounted for by passing Tadj into the inner function, area assignment
		# can happen either there or on the outside. If it happens inside, that info would need 
		# to come from Tadj. If it happens out here the info can come from Pop. For some reason
		# I prefer that. So I'll add an assignment line out here. Assignment will assume that for each
		# year, there can only be one non-NA Area value, and it populates the NAs with that.
		ECpop             <- assignArea(ECpop, Psex)
		Pold              <- Psex[!with(Psex, (Year - Agei - 1) <= omega["Cohmax"] & 
								Agei >= a), ColnamesKeep]
#  dev.new()
#  LexisMap(acast(Pold, Agei~Year, value.var = "Population"),log=FALSE)
		PopMF[[Sex]]      <- rbind(Pold, ECpop[, ColnamesKeep])
	}
	
	if (UNKTF){
		PopMF[["UNK"]]    <- UNK
	}
	
	# stick together males and females, resort and return
	Pop                 <- resortPops(do.call(rbind, PopMF))
	rownames(Pop)       <- NULL
	invisible(Pop)
}

#'
#' @title p_ecm_inner does the work of \code{p_ecm()} 
#' 
#' @description This implementation seeks to follow the MP, we have outsourced the work of \code{p_ecm()} to here for the sake of modularity. This function may also be called by the \code{p_ic()} function ecosystem in order to extend population counts in census 1 and census 2 prior to interpolating. Deaths must be finalized. This function covers area B from Figure 6 'Methods used for population estimates'. 
#' 
#' @param Dsex Deaths after all processing is done, subset of a single sex.
#' @param Tsex Territorial adjustment file, containing only \code{Vx} values for a single sex. Could be \code{NULL}.
#' @param a lowest age for which EC estimates should be given. Default 80.
#' @param omega is the object returned by \code{p_ecm_findOmega()}.
#' @param reproduce.matlab logical. default is \code{FALSE}. This affects only the border cohort between EC and SR.
#' 
#' 
#' @return Population with old ages either filled in or imputed using the extinct cohort method
#' 
#' @details Conceivably, a population could have a long series with low open ages in the early part, but decent population data in the recent part. If we want a value of \code{A} lower than 80 for only part of the series, this can be achieved by subsetting \code{Dsex}, and making two calls to the function, then \code{rbind()}ing back together.
#' 
#' @importFrom reshape2 acast
#' 
#' @export
#' 

p_ecm_inner <- function(Dsex, Tsex = NULL, a = 80, omega = NULL, reproduce.matlab = FALSE){
	
	# instead of feeding in Pop argument, just assume standard R LexisDB format of Pop DF.
	# this is what we do anyway by filling out columns below.
	Pop <- structure(list(PopName = character(0), Area = character(0), Sex = character(0), 
					Age = character(0), AgeInterval = character(0), Type = character(0), 
					Day = integer(0), Month = integer(0), Year = integer(0), 
					RefCode = character(0), Access = character(0), Population = numeric(0), 
					NoteCode1 = character(0), NoteCode2 = character(0), NoteCode3 = character(0), 
					LDB = integer(0), Agei = integer(0), AgeIntervali = integer(0)), .Names = c("PopName", 
					"Area", "Sex", "Age", "AgeInterval", "Type", "Day", "Month", 
					"Year", "RefCode", "Access", "Population", "NoteCode1", "NoteCode2", 
					"NoteCode3", "LDB", "Agei", "AgeIntervali"), row.names = integer(0), class = "data.frame")
	
	if (!"Cohort" %in% colnames(Dsex)){
		Dsex            <- d_addCohortColumn(Dsex)
	}
	if (is.null(omega)){
		omega             <- p_ecm_findOmega(Dsex, l = 5, threshold = 0.5)
	}
	
	
	# 2) reshape deaths by period-cohort (PC = VV shape)
	# VV matrix has years in rows, cohorts in columns
	VV                <- acast(Dsex[with(Dsex, Agei >= a & Cohort <= omega["Cohmax"]), ], 
			Year ~ Cohort, sum, value.var = "Deaths")
	VV                <- rbind(VV, 0)   
	rownames(VV)[nrow(VV)] <- max(Dsex$Year) + 1
	
	# TR: new Aug 25, 2016 for Tadj
	# Assign Vx by cohort of deaths. Problem: each triangle will have one, so each 
	# VV will have 2. Solution acast with function unique().
	if (! is.null(Tsex)){
		# same dims as VV...
		TVV               <- acast(Tsex[with(Tsex, 
								Age >= a & 
								Cohort <= omega["Cohmax"] & 
								Cohort %in% Dsex$Cohort), ],
				  Year ~ Cohort, value.var = "Value", fill = 1)
		# shift just like we shift VV
		TVV              <- rbind(TVV[-1, ],1)
	} else {
		# make 1s if there isn't a Tadj
		TVV               <- VV * 0 + 1
	}
		
	if (reproduce.matlab){
		VV[nrow(VV),ncol(VV)] <- VV[nrow(VV),ncol(VV)] + omega["Ds"]
	}
	cohorts           <- as.integer(colnames(VV))
	years             <- as.integer(rownames(VV))
	
	# all possible ages given these years and cohorts (same dim as VV)
	AllAges           <- outer(years, cohorts, "-") - 1
	# years repeated, same dim as VV, for selection
	Allyears          <- replicate(length(cohorts), years)
	Allcohorts        <- t(replicate(length(years), cohorts))
	# the selection mask
	Mask              <- AllAges >= a & AllAges <= 130
	
	# this may need to iterate

    # TR: new Aug 25, 2016 for Tadj
	# TVV1 conforms with VV, and it will just be 1s for non-tadj instances.
    # this code chunk is always run. No need for more if statements.
	TVVdo             <- apply(TVV[nrow(TVV):1, ], 2, cumprod)[nrow(TVV):1, ]
	TVVundo           <- apply(1 / TVV[nrow(TVV):1, ], 2, cumprod)[nrow(TVV):1, ]
	
	# first scale deaths to same universe
	VVT               <- TVVdo * VV
	# 3) Equation 38
	# Then get pop estimates within single universe
	ECpops            <- apply(VVT[nrow(VVT):1, ], 2, cumsum)[nrow(VVT):1, ]
	# now scale all pops down to their own universe...
	ECpops            <- ECpops * TVVundo
	# use selection mask to just keep ages >= 80, etc.
	ECpops            <- ECpops[Mask]
	
	# remainder of code just to put everything in a standard format.
	# now select the corresponding ages and years
	ECAges            <- AllAges[Mask]
	ECYears           <- Allyears[Mask]
	ECCohorts         <- Allcohorts[Mask]
	
	# stick together and add remaining columns 
	ECpop             <- as.data.frame(
			               matrix(nrow=length(ECAges), 
					              ncol = ncol(Pop), 
					              dimnames = list(NULL, colnames(Pop))))
				  
	ECpop$Age         <- ECpop$Agei            <- ECAges
	ECpop$Year        <- ECYears
	ECpop$Population  <- ECpops
	ECpop$Sex         <- unique(Dsex$Sex)
	ECpop$AgeInterval <- ECpop$AgeIntervali    <- 1
	ECpop$PopName     <- unique(Dsex$PopName)
	
	ECpop$Day         <- 1
	ECpop$Month       <- 1
	ECpop$LDB         <- 1
	ECpop$Cohort      <- ECCohorts
	ECpop$Access      <- "O" # presumably any data we invent is open access, even if the origin data are not?
	ECpop             <- assignNoteCode(ECpop, "p_ecm()")
	ECpop[,colnames(Pop)]
}
