proc import datafile="\yourfilepath"
        out=cox
        dbms=csv
        replace;
    

     getnames=yes;
run;

proc sort data=cox;by simulation N;run;

%macro ODSOff(); 
ods graphics off;
ods exclude all;
ods noresults;
%mend;
 
%macro ODSOn(); 
ods graphics on;
ods exclude none;
ods results;
%mend;

%ODSOff
proc mcmc data=cox seed=801 nbi=1000 nmc=9000 thin=5;
    by simulation N;
    parm mu tau;
    prior mu ~ normal(0, var=16);
    prior tau ~ gamma(shape=.001, iscale=.001);
    random mui ~ normal(mu, prec=tau) subject=cohort monitor=(mui);
    model mean ~ n(mui, var=sd);
	ods output  PostSumInt=post;
run;
%ODSOn

proc export data=post outfile="C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\bhm_unequal_ss_1.csv"
        dbms=csv
        replace;

run;
