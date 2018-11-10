## run matlab script 
# load SWE;
# 
# % d=selif(deaths, deaths(:, end)~=0);  % MISTAKE
# deaths=selif(deaths, deaths(:, end)~=0);
# deaths=deaths(:,1:end-1);
# d1=selif(deaths, deaths(:,10) ~= 21 | deaths(:,3)<1861);
# d2=selif(deaths, deaths(:,10) == 21 & deaths(:,3)>=1861);
# 
# d=d_unk5(d1, d2); 
# save('d1pre.mat', 'd1')
# save('d2pre.mat', 'd2')
# save('dafter.mat', 'd')

library(R.matlab)
library(devtools)
#devtools::install_github("IDButils","timriffe",subdir="IDButils/IDButils")
#devtools::install_github("LDButils","carlboe",subdir="LDButils/LDButils")
#d1mat <- readMat("InputDB/d1pre.mat")
d1mat <- loadMdeaths2R("InputDB/d1pre.mat", "SWE", keepTOT=TRUE)
d2mat <- loadMdeaths2R("InputDB/d2pre.mat", "SWE", keepTOT=TRUE)
daftermat <- loadMdeaths2R("InputDB/dafter.mat", "SWE", keepTOT=TRUE)

isSundbarg <- (Deaths$RefCode ==21 & Deaths$Year >= 1861)  # CHECK! Ref says plenty of contributions below year 1861, so should not matter
d1 <- Deaths[ ! isSundbarg , ]
d2 <- Deaths[   isSundbarg , ]

write.csv(d1, file="d1.csv", row.names=FALSE)
write.csv(f.NC2char(d1mat), file="d1mat.csv", row.names=FALSE)
write.csv(d2, file="d2.csv", row.names=FALSE)
write.csv(f.NC2char(d2mat), file="d2mat.csv", row.names=FALSE)

d.ans=d_unk5(d1, d2);

## compare d.ans with daftermat
d.ans <-d.ans[ order(d.ans$Year, d.ans$Sex, d.ans$Age),]
daftermat <-daftermat[ order(daftermat$Year, daftermat$Sex, daftermat$Age),]
## consider mods to loadMdeaths2R to go against documentation and output NoteCodeX as.character
f.NC2char <- function(x){
  x$NoteCode1 <- as.character(x$NoteCode1)
  x$NoteCode2 <- as.character(x$NoteCode2)
  x$NoteCode3 <- as.character(x$NoteCode3)
  return(x)
}
daftermat <- f.NC2char(daftermat)
head(daftermat)
head(d.ans)
## different # of rows, so look a diff of output
write.csv(daftermat, file="daftermat.csv", row.names=FALSE)
write.csv(d.ans, file="d.ans.csv", row.names=FALSE)


