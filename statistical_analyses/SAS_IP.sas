/*Stool IP study in VLBW babies 2020*/
/*Step 1:Data importing*/
filename stool_IP'/folders/myfolders/SAS Training/stool IP study/cross_ip_dataB.xlsx'; 
options nofmterr;

proc import datafile=stool_IP out=cross_A dbms =xlsx  replace;
run;

/* Managing Variables*/
/* race and sex */
proc format;
	value sex 1 = "Male"
	          2 = "Female";
run;

data cross_A1;
	set cross_A;
	if race = "White" then white_race = 1; else white_race = 0;	
	if race = "AA" then AA_race = 1; else AA_race = 0;	
	if race = "other" or race = "Hispanic" or race = "Asian" then other_race =1; else other_race= 0;
	format sex sex.;
	/*IP_at_1_week = input(IP_at_1_week, best21.);*/
run;

data cross_A2;
	set cross_A1;
	if race = "White" then race2 = "White";	
	if race = "AA" then race2 = "AA";
	if other_race = 1 then race2= "Other";
run;

/* ip*/
proc format;
  	value IP_gut 0 = "Low IP"
  			   	 1 = "High IP";
run;
 			
data cross_A3;
	set cross_A2;	
	if IP_at_1_week > = 0.05 then IP_gut = 1;
	if IP_at_1_week < 0.05 then IP_gut =  0;
	if IP_at_1_week = . then IP_gut = .;
	format IP_gut IP_gut.;
run;

/* birthweight*/
proc format;
  	value bw_cat 0 = "Very Low Birth weight"
  			   	1 =  "Low birth weight";
run;
 
 
data cross_A4;
	set cross_A3;
	if BW <  1500 then bw_cat=0;
	if BW > = 1500 then bw_cat=1;
	if BW  = . then bw_cat=. ;
	format bw_cat bw_cat.;
	
	if PMA < 31 then PMA_cat = 0;
	else PMA_cat=1;
run;

/*gestational age*/

proc format ;
	value GA_cat 0 = "GA ≤ 28weeks"
			     1 = "GA > 28weeks";
run;

data cross_A5;
	 set cross_A4;
	 if GA < = 28 then GA_cat = 0;
	 if GA >  28 then GA_cat = 1;
	 if GA < = . then GA_cat = .;
	 format GA_cat GA_cat.;
run;
/* baby antibiotic days*/

proc format;
	value baby_3abdays 0 = "Received abx ≤ 3 days"
					   1 = "Received abx > 3 days";		    
run;	 

proc format; 
	value delroute 1 ="Vaginal"
	               2 = "Cesarean";
run;
	
data cross_A6;
	set cross_A5;
	
	if baby_abdays < = 3 then baby_3abdays = 0;
    if baby_abdays >  3 then baby_3abdays = 1;	
	
	format baby_3abdays baby_3abdays.;
	format delroute delroute.;
run;
		 
/* Mothers breast milk days*/

proc format;
	value MBM_3Feed 1= "MBM for < 3 days before IP test"
                    0= "MBM for > = 3 days before IP test" ; 
run;

data cross_A7;
	set cross_A6;
    
    if numDayMBMFeed < 3 then MBM_3Feed = 0 ;
    if numDayMBMFeed > = 3 then MBM_3Feed = 1;
    
    format MBM_3Feed MBM_3Feed.;
   
run;    

/* Combined mothers breast milk feed days*/

proc format;
	value varF 0 = "Other combinations"
			   1 = "Antibiotics < = 3 and MBM > = 3";
run;

data cross_A8;
 	set cross_A7; 
  
    if baby_abdays < = 3 and numDayMBMFeed > = 3 then varF = 1;
    ELSE varF = 0;
    format varF varF.;
    
run;

proc contents data= cross_A8 position;
run;

proc export data=cross_A8 outfile= '/folders/myfolders/SAS Training/stool IP study/final_cross_exported.xlsx'
	dbms= xlsx replace;
run;

/*Descriptive statistics */

proc means data=cross_A8 n mean stddev median min max maxdec=2;
 	var BW  GA  PMA  Apgar1 Apgar5  ;
run;

proc freq data=cross_A8;
	tables sex race2 delroute PPROM pre_eclampisa antsteroids IP_gut bw_cat GA_cat 
	mat_abs baby_abdays numDayMBMFeed baby_3abdays MBM_3Feed varF PMA_cat;			
run;
	 
/*table 1A stratified by IP using fisher's exact test*/

data cross_IP_chi;
	set cross_A8;
run;

proc freq data=cross_IP_chi;
	 tables(sex race2 white_race AA_race other_race delroute PPROM pre_eclampisa antsteroids mat_abs baby_abdays bw_cat
	 GA_cat baby_3abdays MBM_3Feed varF PMA_cat)* IP_gut/ fisher nocol;
run;
	
	 
/* table 1A continous variables*/

data cross_ttest;
	set cross_A8;
run;

proc ttest data = cross_ttest;
	class IP_gut;
	var GA;
run;	

proc ttest data = cross_ttest;
	class IP_gut;
	var Apgar1;
run;

proc ttest data = cross_ttest;
	class IP_gut;
	var Apgar5 ;
run;


proc ttest data = cross_ttest;
	class IP_gut;
	var  PMA ;
run;

proc ttest data = cross_ttest;
	class IP_gut;
	var BW  ;
run;

/*Table 1B unadjusted  and adjusted odds ratio for factors associated with IP*/
data cross_logistic;
	set cross_A8;
run;

proc logistic data=cross_logistic;
  	class baby_3abdays(ref="Received abx > 3 days") / param=glm ;
  	model IP_gut( event= "High IP")= baby_3abdays ;
run;

proc logistic data=cross_logistic;/* adjusted*/
  	class baby_3abdays(ref="Received abx > 3 days") / param=glm ;
  	model IP_gut( event= "High IP")= baby_3abdays PMA BW;
run;

proc logistic data=cross_logistic;
  	class MBM_3Feed (ref="MBM for < 3 days before IP test") / param=glm ;
  	model IP_gut( event= "High IP")= MBM_3Feed  ;
run;

proc logistic data=cross_logistic;/* adjusted*/
  	class MBM_3Feed (ref="MBM for < 3 days before IP test") / param=glm ;
  	model IP_gut( event= "High IP")= MBM_3Feed BW PMA  ;
run;

proc logistic data=cross_logistic;
  	class varF (ref= "Other combinations") / param=glm ;
  	model IP_gut( event= "High IP")= varF ;
run;

proc logistic data=cross_logistic;/* adjusted*/
  	class varF (ref= "Other combinations") / param=glm ;
  	model IP_gut( event= "High IP")= varF BW PMA;
run;

/* To check for collinearity */

proc corr data=cross_A8 spearman;
var PMA numDayMBMFeed baby_abdays BW GA sex delroute
PPROM pre_eclampisa antsteroids mat_abs Apgar1 Apgar5
white_race AA_race other_race IP_gut VarF MBM_3Feed baby_3abdays ;
title 'IP study variables - Examination of Correlation
Matrix';
run;

	 
proc reg data=cross_A8;
model IP_gut = PMA numDayMBMFeed baby_abdays BW GA sex delroute
PPROM pre_eclampisa antsteroids mat_abs Apgar1 Apgar5
white_race AA_race other_race/ vif tol collin;
title 'Multicollinearity Investigation of VIF and Tol';
run;
quit;


