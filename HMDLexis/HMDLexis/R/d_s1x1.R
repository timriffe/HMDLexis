
#'
#' @title d_s1x1, a function to split RR death counts into Lexis triangles based on regression coefficients.
#' 
#' @description The HMD Methods Protocol gives a formula to split Lexis 1x1 death counts (squares) into triangles based on the results of a regression. Function can be run innocuously on any Deaths data, even if no 1x1 RR is present.
#' 
#' @details This function can optionally deal with territorial adjustments. If a Tadj file is given, it is handled appropriately. Since Tadj files are created by default by \code{readInputDB()}, there is no reason not to specify them always, even when not relevant. It doesn't matter. If you have a Tadj, though, you better know it.
#' 
#' @param Deaths LexisDB internal Deaths \code{data.frame}, long format, all columns. Format as given by \code{readInputDB()}.
#' @param Births LexisDB internal Births \code{data.frame}, as given by \code{readInputDB()}.
#' @param Tadj LexisDB internal territorial (or universe) adjustment \code{data.frame}, as given by \code{readInputDB()}. This is optional.
#' 
#' @importFrom reshape2 acast
#' @importFrom reshape2 melt
#' @importFrom compiler cmpfun
#' 
#' @export 
#' 

d_s1x1 <- function(Deaths, Births, Tadj = NULL){
	# check if function needs to be run, if not, return deaths 
	# TR: 1 June, 2016: use new %==%, more efficient
	if (!any(Deaths$Lexis[with(Deaths, AgeIntervali %==% 1 & 
									YearInterval == 1)] == "RR")){
		cat("d_s1x1() not necessary; no 1x1 RR deaths to split into triangles at this time.")
		return(Deaths)
	}  
	# some prelim data procedures, pretty run-of-the-mill
	#--------------------------------------------------
	# TOT is never necessary:
	Deaths <- Deaths[Deaths$Age != "TOT", ]
	
	# slice off UNK, rbind back on later:
	UNKi <- Deaths$Age == "UNK"
	UNKTF <- any(UNKi)
	if (UNKTF){
		UNK           <- Deaths[UNKi, ]
		Deaths        <- Deaths[!UNKi, ]
	}
	
	# slice off OP, to be slapped back on later
	OPi <- Deaths$AgeInterval == "+" 
	OPTF <- any(OPi)
	if (OPTF){
		OP            <- Deaths[OPi, ]
		Deaths        <- Deaths[!OPi, ]
	}
	
	# TR: added Aug 22, 2016
	TadjTF <- is.null(Tadj)
	#--------------------------------------------------
	# now start the sex-loop (this can't be done in a big Year-Sex apply scheme
	# because consecutive birth cohort sizes are needed
	Areas <- unlist(tapply(Deaths$Area, Deaths$Year,unique))
	Dout  <- list()
	# Sex <- "f"
	for (Sex in c("f","m")){
		Dsex            <- Deaths[Deaths$Sex == Sex, ]                                                  
		RRi             <- with(Dsex, Lexis == "RR" & AgeInterval == "1" & YearInterval == 1)
		RR              <- Dsex[RRi, ]
		DO              <- Dsex[!RRi, ] # Save other Lexis shapes for downstream treatment.
		
		maxAges         <- unlist(lapply(split(RR, RR$Year), function(RRyr){
							max(RRyr$Agei)
						}))
		
		# note that this does the job of d_ma0() for a large block of ages
		RRwide          <- acast(RR, Agei ~ Year, sum, value.var = "Deaths", fill = 0)
		# TR: Aug 23, 2016. Note that RRwide needn't be a continuous block of years.
		# years can have gaps, since it only picks up years with RR to split...
		
		# TR: 1 June, 2016
		# ensure that RRwide goes to age 0
		AgesBox         <- 0:as.integer(rownames(RRwide)[nrow(RRwide)])
		yrs             <- as.integer(colnames(RRwide))
		RRwideBox       <- matrix(0,
				ncol=length(yrs),
				nrow=length(AgesBox),
				dimnames=list(AgesBox, yrs))
		RRwideBox[rownames(RRwide),colnames(RRwide)] <- RRwide
		RRwide          <- RRwideBox
		# years are in order, but might not be continuous chunks..........
		yrsc 			<- colnames(RRwide)
		yrs             <- as.integer(yrsc)
		
		# TR: modified Mon Aug 22, 2016 to account for Tadj.
		# begin modify here.
		Ball            <- Births$Births[Births$Sex == Sex]
		names(Ball)     <- Births$Year[Births$Sex == Sex]
		Ball            <- Ball[sort(names(Ball))]
		
		BT_1 			<- Ball[as.character(yrs - 1)]
		BT   			<- Ball[yrsc]
		if (!TadjTF){
			ind 		<- Tadj$Type == "Rb" & Tadj$Sex == Sex
			RB 			<- Tadj$Value[ind]
			names(RB) 	<- Tadj$Year[ind]
			# need names BT...
			BT_1        <- BT_1 * RB[names(BT)]
		}
		IMRdenom        <- (1 / 3 * BT_1) + (2 / 3 * BT)
		# names can be shifted depending on order of BT_1 and BT. Might be relevant for alignment
		names(IMRdenom) <- names(BT)
		# end modify here
		# this would only kick in if births from year t-1 are missing. Then assume same as present year
		NAind           <- is.na(IMRdenom)
		IMRdenom[NAind] <- Ball[yrsc[NAind]]
		
		# this is trickier than first glance, since it's conceivable to have infant mort
		# in TL TU while higher ages are in RR, which need to be split. We therefore don't
		# necessarily take the first row of RRwide as infant mort, but rather sum it independently
		# VH, VV, or RV mort is ignored here.
		
		# TR: Aug 22, 2016. Noted, numerator does not need Tadj adjust, only denom.
		IMRnum          <- acast(Dsex[Dsex$Lexis %in% c("TL","TU","RR") & Dsex$Agei == 0, ], 
				Age ~ Year, 
				sum, 
				value.var = "Deaths", 
				fill = 0)[ , , drop = TRUE] # drop = TRUE turns it into a vector (named still)
		
		IMR             <- IMRnum[yrsc] / IMRdenom # ensures we get same years as need for RR
		
		# get these in a matrix conformable with formulas
		IMRT            <- IMR[col(RRwide)]
		dim(IMRT)       <- dim(RRwide)
		
		# now robust to non-consecutive yrs
		# account for fluctuations in cohort size (needs tadjification)
		# this takes ALL available cohorts, not just those constrained by RR needs
		
# TR: changed Mon Aug 22, 2016 to use Tadj where necessary
# begin modify here
		BT_1 			<- Ball[-length(Ball)]
		BT   			<- Ball[-1]
		if (!TadjTF){
			# similar lines, but possibly different yrIn
			#yrIn  <- as.character(as.integer(names(BT_1)) + 1)
			BT_1  <- BT_1 * RB[names(BT)]
		}
		pib             <- BT / (BT + BT_1)
#   pib             <- Ball[2:length(Ball)] / 
#                          (Ball[2:length(Ball)] + Ball[1:(length(Ball) - 1)])
		## end modify here                                 
		# now we determine the cohort for each cell in RRwide :-)
		yrsM            <- yrs[col(RRwide)]
		ageM            <- as.integer(rownames(RRwide))[row(RRwide)]
		TLCoh           <- as.character(yrsM - ageM)
		dim(TLCoh)      <- dim(RRwide)
		# now distribute these over cohorts
		PIB             <- pib[TLCoh]
		dim(PIB)        <- dim(RRwide)
		dimnames(PIB)   <- dimnames(RRwide)
		# assume constant sizes for cohorts not available
		# plot(density(PIB, na.rm=TRUE))
		PIB[is.na(PIB)] <- 0.5
		# [[ NOTE: if births prior to first year of deaths are available, we could use these too ]]
		
		# some indicator matrices
		Ix0 <- Ix1 <- I1919 <- I1918 <- PIB * 0
		if (1918 %in% yrs){
			I1918[,"1918"] <- 1
		}
		if (1919 %in% yrs){
			I1919[,"1919"] <- 1
		}
		Ix0["0", ]       <- 1
		Ix1["1", ]       <- 1
		
		if (Sex == "f"){
			# from Table A1, ages 0-130
			alpha <- c(0.0392, 0.1365, rep(0.0130,3),rep(c(0.0018,-0.0140,-0.0135,-0.0061,-0.0046,-0.0041,-0.0072,-0.0070,
									-0.0071,-0.0084,-0.0091,-0.0134,-0.0175,-0.0201,-0.0230,-0.0231,-0.0187,-0.0112,-0.0014), each = 5),
					rep(0.0190,31))
			names(alpha) <- 0:130
			alpha <- alpha[rownames(RRwide)]
			ALPHA <- (PIB * 0 + 1) * alpha
			
			# the full formula
			pTL <-  0.4710 + ALPHA + 0.7372 * (PIB - 0.5) +
					0.1025 * I1918 -0.0237 * I1919 +
					-0.0112 * log(IMRT) +
					-0.0688 * log(IMRT) * Ix0 +
					0.0268 * log(IMRT) * Ix1 +
					0.1526 * (log(IMRT) - log(0.01)) * Ix0 * (IMRT < 0.01)
		}
		if (Sex == "m"){
			# from Table A1, ages 0-130
			alpha <- c(0.0230, 0.1249, rep(0.0086, 3), rep(c(0.0031, -0.0086, -0.0175, 0.0035, 0.0081, 0.0031,
									-0.0065, -0.0117, -0.0148, -0.0145, -0.0142, -0.0157, -0.0179, -0.0198, -0.0223, -0.0216,
									-0.0160, -0.0083, 0.0039), each = 5), rep(0.0313, 31))
			names(alpha)  <- 0:130
			alpha         <- alpha[rownames(RRwide)]
			ALPHA         <- (PIB * 0 + 1) * alpha
			# the full formula
			pTL           <-  0.4838 + ALPHA + 0.6992 * (PIB - 0.5) +
					0.0728 * I1918 -0.0352 * I1919 +
					-0.0088 * log(IMRT) +
					-0.0745 * log(IMRT) * Ix0 +
					0.0259 * log(IMRT) * Ix1 +
					0.1673 * (log(IMRT) - log(0.01)) * Ix0 * (IMRT < 0.01)
		}
		# get counts in triangles
		TL              <- RRwide * pTL
		TU              <- RRwide - TL
		
		# put in long format
		TL              <- melt(TL, varnames = c("Age","Year"), value.name = "Deaths")
		TU              <- melt(TU, varnames = c("Age","Year"), value.name = "Deaths")
		
		# select only ages up until original max age of RR in given year
		TL              <-  do.call(rbind,
				lapply(
						split(TL, TL$Year), function(YR, maxAges){
							yr  <- unique(YR$Year)
							YR[YR$Age <= maxAges[as.character(yr)],]
						},maxAges=maxAges)
		)
		TU              <- do.call(rbind,
				lapply(
						split(TU, TU$Year), function(YR, maxAges){
							yr <- unique(YR$Year)
							YR[YR$Age <= maxAges[as.character(yr)],]
						},maxAges=maxAges)
		)
		
		# add Lexis shape
		TL$Lexis           <- "TL"
		TU$Lexis           <- "TU"
		
		DeathsTLTU         <- rbind(TL, TU)
		
		DN                 <- as.data.frame(
				matrix(ncol = ncol(Deaths), 
						nrow = nrow(DeathsTLTU), 
						dimnames = list(NULL, colnames(Deaths)))
		)
		
		DN$PopName         <- unique(Deaths$PopName)
		DN$Year            <- as.integer(DeathsTLTU$Year)
		DN$Area            <- Areas[as.character(DN$Year)]
		DN$YearInterval    <- 1
		DN$Sex             <- Sex
		DN$Age             <- as.character(DeathsTLTU$Age)
		DN$AgeInterval     <- "1"
		DN$Lexis           <- DeathsTLTU$Lexis
		DN$Deaths          <- DeathsTLTU$Deaths
		DN$Agei            <- DeathsTLTU$Age
		DN$AgeIntervali    <- 1
		DN$LDB             <- 1
		DN                 <- assignNoteCode(DN, "d_s1x1()")
		DN$Access          <- "O" # this maybe need to be sensitive? apply general rule later?
		Dout[[Sex]]        <- rbind(DO, DN)
	}
	
	# tack on parts we sliced off
	if (UNKTF){
		Dout[["UNK"]]      <- UNK
	}
	if (OPTF){
		Dout[["OP"]]       <- OP
	}
	Deaths               <- resortDeaths(do.call(rbind, Dout))
	
	# TR: 1 June, 2016: add d_agg() step in case RR overlapped with TL,TU
	Deaths               <- d_agg(Deaths)
	rownames(Deaths)     <- NULL
	invisible(Deaths)
}



