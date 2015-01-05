#'
#' @title d_mao a function to add missing ages in deaths with zero counts
#' @description The function loops a lot to fill missing age-lexis combos with 0s
#' 
#' @return Deaths the standard Deaths data.frame, with 0s imputed as necessary

#Lexis element (2-char). This field denotes the shape of the Lexis
#element, where: TL=lower triangle, TU=upper triangle, RR=rectangle,
#VV=parallelogram with vertical left and right sides (i.e.,
#period-cohort), VH=parallelogram with horizontal upper and lower
#sides (i.e., age-cohort); RV= Same as VV except also includes TL for
#the first age in the interval (e.g., cohorts aged 0-4 on Dec 31st
#—includes those born in the current calendar year

## the entirety of Lexis shapes appears in the collection (DNK,FIN,BEL). The RV and VH shapes
## are rare and are handled by d_ma01()

d_ma0 <- function(Deaths, sa=90){
  ## adds missing ages in deaths with deaths numbering 0
  ## Input: deaths array
  ## Output: deaths array
  
  ## for each year, sex, see what ages are covered by current ages and intervals,
  ## determine if rows missing assuming 0:omega as ages, fill in.
  
  ## This could probably be made much faster and easier to debug by (1) separating the sexes (and remove the condition for change of sex across rows); (2) chunking by year -  there is then no need to 
  
  
  d<-Deaths
  d<- d[order(d$Sex, d$Year, d$Age, d$Lexis),]
  
  daccum <- d[0,];
  
  
  for( sex in unique(d$Sex) ){
    isel.sex <- (d$Sex == sex)
    for( yr in unique(d$Year) ){
      #print(paste("sex, year", sex, yr))
      
      isel.year <- (d$Year == yr)
      
      d.sexyr <- d[(isel.sex & isel.year),] # only within sex, year
      
      ## split the data into two pieces, young|old
      isel.age <- d.sexyr$Agei < sa | d.sexyr$Agei == 300
      d1 <- d.sexyr[(isel.age),]
      d2 <- d.sexyr[!(isel.age),]
      ## recode $Lexis, " 1='TL'; 2='TU'; 3='RR'; 6='RV'; 4='VV' ; 5='VH'; -1=NA "
      
      for( i in seq_len(nrow(d2))[-nrow(d2) ] ){
        if (d2$AgeIntervali[i] > 1 ){
          d1 <- rbind(d1, d2[i,] )
          next
        }
        
        ## vars to make readable
        t.iVV <- d2$Lexis[i] =='VV'
        t.ip1VV <- d2$Lexis[i+1] =='VV'    
        t.iTL <- d2$Lexis[i] =='TL'
        t.ip1TU <- d2$Lexis[i+1] =='TU'
        
        t.ageeq <- d2$Agei[i+1] == d2$Agei[i]
        t.deltaage1 <- d2$Agei[i+1] == d2$Agei[i] + 1
        t.yreq  <- d2$Year[i] == d2$Year[i+1]
        t.sexeq <- d2$Sex[i] == d2$Sex[i+1]
        
        if(t.iVV & t.ip1VV  & !t.deltaage1 ){ # check this!
          ## this is selected if VV data and i is not the final VV in the year.
          ## adjacent VV with AgeInterval more than 1 i.e. d2$AgeIntervali > 1
          d1 <- rbind(d1, d2[i,])
          d2$Deaths[i] <- 0
          d2$Agei[i] <- d2$Agei[i] + 1 ;
          d2$Age[i] <- as.character( d2$Agei[i]) # parallel char and num Age
        } else if(  t.iTL  & ( (t.ip1TU  & t.ip1VV & t.ageeq  & t.yreq ) |  d2$Agei[i] + 1 < d2$Agei[i+1] )  ) {
          ## potential TL --> TU transformation
          d1 <- rbind(d1, d2[i,])
          d2$Deaths[i] <- 0
          d2$Lexis[i] <- 'TU'
        }
        
        if( d2$Agei[i] == 1 & d2$Agei[i+1] == 2 &  t.sexeq & t.yreq  & t.ageeq ){
          d1 <- rbind(d1, d2[i,])
          d2$Deaths[i] <- 0
          d2$Lexis <- 'TU'
        }
        
        ##%  if d(i,4)+1<d(i+1,4) & d(i,2)==d(i+1,2) & d(i,3)==d(i+1,3) 
        if( (d2$Agei[i] + 1 < d2$Agei[i])  &  t.sexeq & t.yreq ){
          ##for k=1:d(i+1,4)-d(i,4)-1
          for(k in 1:(d2$Agei[i+1] - d2$Agei[i] - 1) ){ # watch out! matlab vs R precedence on ':'
            t.Nd1 <- nrow(d1)               #length(d1(:,1)
            d1 <- rbind(d1, d2[i,]);       
            d1$Agei[nrow(d1)] <- d$Agei[i]+k
            d1$Age[nrow(d1)] <- as.character(d1$Age[nrow(d1)])
            d1$Deaths[nrow(d1)] <- 0;          #d1(length(d1(:,1)),9)=0;
            d1$Lexis[nrow(d1)] <- 'TL';        #d1(length(d1(:,1)),5)=1;
            d1<- rbind(d1,d2[i,]);          #d1=[d1; d(i,:)]; 
            d1$Agei[nrow(d1)]<- d2$Agei[i]+k; #d1(length(d1(:,1)),4)=d(i,4)+k;
            d1$Age[nrow(d1)] <- as.character(d1$Agei[nrow(d1)])
            d1$Deaths[nrow(d1)]<- 0;          #d1(length(d1(:,1)),9)=0;
            d1$Lexis[nrow(d1)]<-'TU';         #d1(length(d1(:,1)),5)=2;
          }
          
          if( d2$Lexis[i] =='TU' | d2$Lexis[i] == 'VV') { #if d(i+1,5)==2 | d(i+1,5)==4
            d1 <- rbind(d1, d2[i,]);        #d1=[d1; d(i+1,:)];
            d1$Deaths[nrow(d1)]<- 0;           #d1(length(d1(:,1)),9)=0;
            d1$Lexis[nrow(d1)] <- 'TL';        #d1(length(d1(:,1)),5)=1;
          }
          
        }
        
        d1 <- rbind(d1, d2[i,]);            #d1=[d1; d(i,:)];
        
        
        if (d2$Lexis[i]=='TU' & d2$Lexis[i+1] == 'TU' & t.ageeq){ #if d(i,5)==2  & d(i+1,5)==2 & d(i,4)+1==d(i+1,4)
          d1 <- rbind(d1, d2[i+1,]);        #d1=[d1; d(i+1,:)]
          d1$Deaths[nrow(d1)] <- 0;       #d1(length(d1(:,1)),9)=0;            
          d1$Lexis[nrow(d1)]<- 'TL';       #d1(length(d1(:,1)),5)=1; 
        } 
        
      }
      
      d1 <- rbind(d1, d2[nrow(d2),]);        #  d1=[d1; d(length(d(:,1)),:)];
      
      daccum <- rbind(daccum, d1);
      #print(paste("size of daccum:",nrow(daccum)) )
      
    }  # yr
    
    
  } #sex
  
  
  return(daccum)
}