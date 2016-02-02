# local - country specific defs, included in main Makefile

COUNTRY=HUN
countrylc=hun
CountryName=Hungary


XCOUNTRY=$(COUNTRY)
# enable the following line to allow recomputation of an XYZ_Old/
# collection or for recomputation to occur in a directory other than
# $(PROJ)/$(COUNTRY)
# XCOUNTRY=HUN_Old

#who to notify of Working changes; Country Specialist at the minimum
#Dana Glei (danaglei@pacbell.net)

CS=Jasilionis@demogr.mpg.de,boe@demog.berkeley.edu,Jdanov@demogr.mpg.de, lyang@demog.berkeley.edu, magali@demog.berkeley.edu

# who to notify for publishing; CS,carl,dana,john, vladimir at minimum
NotifyList=$(CS),danaglei@demog.berkeley.edu,jrw@demog.berkeley.edu,Shkolnikov@demogr.mpg.de, lyang@demog.berkeley.edu, magali@demog.berkeley.edu

