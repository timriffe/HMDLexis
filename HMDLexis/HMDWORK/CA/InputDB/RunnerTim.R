
setwd("/data/commons/cwinant/localDataTest")
## setwd("/data/commons/triffe/LexisDB_R/testing/TestStates/")
devtools::load_all("/data/commons/triffe/git/HMDLexis/HMDLexis/HMDLexis")

IDB <- readInputDB(WORKING = "CA", 
  XXX = NULL, 
  log.file = "", 
  InputDB.name = "InputDB",
  save.bin = FALSE, verbose = TRUE, strict = FALSE)
Births <- IDB$Births
Deaths <- IDB$Deaths
Pop    <- IDB$Pop
Pop    <- Pop[Pop$Age != "TOT", ]

# Process deaths, easy for USA in general:
Deaths <- d_unk(Deaths)
Deaths <- d_long(Deaths)
Deaths <- Deaths[Deaths$Agei <= 130, ] # age reporting is crayzay

# process population counts:

# 1) Separate 5-year intercensal data for 1970s...
P5_1970s <-  Pop[Pop$RefCode == "60", ] # save these, use them as a standard.

# we can do a lot of preparing now and then go back to rescale as necessary.
# with the caveat that the 1970s 5-year data is mid-year...

# this is the main Pop object we'll prepare
Pop      <- Pop[Pop$RefCode != "60", ]

# redistribute UNK (1970s has no UNK, so no need there)
Pop      <- p_unk(Pop)

# split partial 5-year age groups for intercensal endpoints...
# (the slowest function in the LexisDB, sorry)
Pop      <- p_split(Pop, Deaths, Births)
# FYI this would break if the 1970s 5-year intercensals weren't separated.

# split open age group. This is only needed for 1980 data prior to 
# running p_ic(). For all other years it will also redistribute
# open ages, but this is innocuous for downstream calcs.
Pop      <- p_soai(Pop, Deaths)

# now do intercensals as necessary (1970s, our first draft, so to speak)
# we will rescale these in a moment to match the midyear 5-year age groups.
Pop      <- p_ic(Pop, Deaths, Births)
# recall, the reason why we do an HMD first draft of 1970s intercensals is to get a better distribution
# for single ages, since, ic() comes from both left and right...

# we've now done intercensals for the 1960s and 1970s. The 1960s intercensals are taken as final

# this gets us Jan 1, 1959 and 1960, removes 1960 census.
Pop      <- p_precensal(Pop, Deaths, Births)

# TR: postcensals written quickly, and may need further testing.
# Matlab compare wasn't carried out. Mostly only need to check
# the infant cohort, since it's the most likely place they would
# have gotten lazy in matlab or where I may have messed up
Pop      <- p_postcensal(Pop, Deaths, Births) 
# Pop goes out to Jan 1, 2014 now
################################################################################
# PART II, prepare data that we use as standard for rescaling: official        #
# intercensals in 5-year age groups, but MIDYEAR. Some dance moves required    #
################################################################################
# steps to prep data to scale to
# * we have abridged counts at mid year. We want Jan 1 counts. Don't feel comfortable moving
# 5-year age groups to Jan 1 using p_movedata() because it's sloppy- a simple average? yuck.

# 1) split age groups, still at mid year. (rbind in in 1980s to use as distribution).
# 2) move data to Jan 1 using overly simple but not sinful assumptions, p_movedata()
# 3) then regroups to 5-year age groups for rescaling, using p_agg(), but this is hidden,
#   I put it inside of p_rescale() out of sight, because we don't want that function available.
# 4) then p_rescale().  

# the point for all this dancing is that that in order to split 5-year official intercensals we'd
# prefer to use distributions both on left and right rather than on right only...
# these we get for free using p_ic()... so it's a parlor trick, really.
################################################################################

# do rescale first

# slow. Let's split 5-year age groups before moving to Jan 1
P1_1970s1      <- p_split(rbind(P5_1970s, Pop[Pop$Year == 1980, ]), Deaths, Births)

# then cut off 1980, because 1980 is Jan 1 and 1979 is July 1.
# we only need the Jan 1 estimates for the 1970s, and Jan 1 1979 is boxed in 
# by July 1978 and July 1979.
P1_1970s1      <- P1_1970s1[P1_1970s1$Year < 1980, ]
# now move these data to Jan 1

# this will give Jan 1, 1971 to Jan 1, 1979.
P_1970s        <- p_movedata(P1_1970s1)

################################################################################
# now rescale our HMD intercensals for 1970s
Pop            <- p_rescale(Pop, P_1970s, N = 5)

# finally p_srecm() closes out the population
Pop            <- p_srecm(Pop, Deaths, A = 85)

# check we have ages 0-130 in each year / sex
unique(table(Pop$Agei)) # 112
################################################################################
# Carl can fill this in.
# writeLDB(...)

