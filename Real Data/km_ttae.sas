libname az 'C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\Real Data\Data\azdata';
libname brii 'C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\Real Data\Data\briidata';
libname lilly 'C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\Real Data\Data\lillydata';
libname mp 'C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\Real Data\Data\mpdata';
libname pf 'C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\Real Data\Data\pfdata';
libname vir 'C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\Real Data\Data\virdata';

proc format;
  value grp
    1 = "Treatment"
    2 = "Placebo"
    ;
run;

/* read data */
%macro readdata(in,out);
data &out;
 set &in;
run;
%mend;
%readdata(in=az.ticopublicaz,out=az);
%readdata(in=brii.ticopublicbrii,out=brii);
%readdata(in=lilly.ticopubliclilly,out=lilly);
%readdata(in=mp.ticopublicmp,out=mp);
%readdata(in=pf.ticopublicpf,out=pf);
%readdata(in=vir.ticopublicvir,out=vir);


/* KM plot - time to recovery*/
%macro kmrecov(in,out);
data adtte2;
  set &in;
  if safetyd90 = 1 then cnsr= 0;
  else cnsr= 1;
  aval=t2safetyd90;
  grpn=rangp;
run;

proc sql noprint;
   select max(AVAL) into: maxsur1 separated by "" from ADTTE2;
quit;	

%put maxsur1 = &maxsur1.;
%let maxsur = %sysfunc(ceil(&maxsur1.));
%put maxsur = &maxsur.;

/*---------------Start Statistical analysis---------------------*/

/*Summary+median PFS/CI+unstrata pvalue + survial rate*/ 

ods graphics on;
ods listing close;
ods output CensoredSummary=Summary 
           productlimitestimates=surv
		   homtests=Homtests_NSTRA
		   Survivalplot=atrisk
           Quartiles=Quartiles;

proc lifetest data=ADTTE2 method=km plots=survival(atrisk=0 to &maxsur. by 5)
  timelist=0 to &maxsur. by 5 atrisk ;
  time AVAL*cnsr(1);
  strata /group=grpn TEST=logrank; 
  SURVIVAL OUT=PPLOT STDERR;
run;
ods output close;
ods listing; 

data dummy;
do grpn = 1,2;
	   output;
	end;

 proc sort; by grpn;
run;

proc sql noprint;
  create table tot as
  select count(*) as tot, grpn
  from ADTTE2
  group by grpn
  order by grpn
  ;
quit;

/*summary*/

proc sort data = summary; by grpn; run;

data summary2;
  merge summary(in=a) dummy(in=b);
  by grpn;

  if missing(Total) then Total = 0;
  if missing(Failed) then Failed = 0;

  proc sort; by grpn;
run;

data summary3;
  merge summary2(in=a) tot(in=b);
  by grpn;

  if missing(tot) then delete;

run;

data _null_;
  set summary3;
  length Event $200.;
  if (Failed ^= 0 and tot ^= 0 and cmiss(failed, tot) = 0) then Event=cats(Failed)||' ('||strip(put(100*Failed/tot,5.1))||'%)';
  else if (Failed = 0 and tot ^= 0) then Event= "0"||' ('||strip(put(100*Failed/tot,5.1))||'%)';
  else if (Failed = 0 and tot = 0) then Event= "0";

  if (Failed ^= 0 and tot ^= 0) then eventN = cats(Failed)||"/"||cats(tot);
  else if (Failed = 0 and tot ^= 0) then eventN = cats(Failed)||"/"||cats(tot);
  else if (Failed = 0 and tot = 0) then eventN= "0";

  call symputx('NPT'||cats(grpn),cats(tot));
  call symputx('Event'||cats(grpn),cats(Event));
  call symputx('EventN'||cats(grpn),cats(eventN));

run;

  %put npt1 = &npt1 npt2 = &npt2 event1 = &event1 event2 = &event2;
  %put eventN1 = &eventN1 eventN2 = &eventN2;



/*median PFS + CI*/

proc sort data = Quartiles(where=(percent=50)); by grpn; run;

data Quartiles2;
  merge Quartiles(in=a) dummy(in=b);
  by grpn;

run;

data _null_;
  length Estimate1 LowerLimit1 UpperLimit1 result $200.;
  set Quartiles2;
  if Estimate=.   then Estimate1='NE';
  else Estimate1=cats(put(round(Estimate,0.01),10.2));
  if LowerLimit=. then LowerLimit1='NE';
  else LowerLimit1=cats(put(round(LowerLimit,0.01),10.2));
  if UpperLimit=. then UpperLimit1='NE';
  else UpperLimit1=cats(put(round(UpperLimit,0.01),10.2));
  result=cats(Estimate1)||' ('||cats(LowerLimit1)||', '||cats(UpperLimit1)||')';
  call symputx('MEDCI'||cats(grpn),result);
run;

  %put medci1 = &medci1 medci2 = &medci2;


/*unstrata P value*/

data _null_;
  length pvalue $200.;
  set Homtests_NSTRA;
  if      ProbChiSq=.        then pvalue ='NE';
  else if .<ProbChiSq<0.0001 then pvalue ='<0.0001';
  else if ProbChiSq>0.9999   then pvalue ='>0.9999';
  else                            pvalue =cats(put(round(ProbChiSq,0.0001),11.4));
  call symputx('NSFP',cats(pvalue));
run;

%put pvalue = &nsfp;

/*at risk number*/

data atrisk2;
  set atrisk;

/*  grpn = input(Stratum,best.);*/

  if missing(tAtRisk) then timepoint = ceil(Time);
  if ^missing(tAtRisk) then timepoint = tAtRisk;

  if survival ^= . then survival=survival*100;
  if censored^=. then censored=censored*100;

  if missing(event) then event = 0;

  if cmiss(event,AtRisk) = 0 then atrisk2 = strip(cats(AtRisk))||' ('||strip(cats(event))||')';

  atrisk1 = strip(put(AtRisk,best.));

  rename tatrisk = arxvalue;

  proc sort; by grpn timepoint;

run;

/*survial rate*/

data survial(drop=temp);
  set atrisk2;
  by grpn timepoint;

  retain temp;

  if first.grpn then call missing(temp);
  if Survival^=. then temp = Survival;
  else if Survival=.  then Survival = temp;

  if censored^=. then CNSR_Survival = Survival;


  proc sort; by grpn timepoint;

run;

/*-------draw figure--------*/

proc sql noprint;
  select unique arxvalue into: TimepointN separated by ' ' from survial
  where arxvalue ne .;
quit;

%put **&TimepointN.**; 

/********************output ****************************/

proc template;
  define statgraph plot;
  begingraph/backgroundcolor=white designheight=4.10in designwidth=9.00in /*DataColors=(blue red) DataSymbols=(CIRCLE X)*/ DATACONTRASTCOLORS=(blue red green purple);

  layout lattice/columns=1 rows=3 rowweights=(0.8 0.05 0.15) border=false opaque=true backgroundcolor=white columndatarange=union;

  layout overlay/cycleattrs=true /*walldisplay=(fill)*/
  yaxisopts=(griddisplay=off label="Probability of AE-Free (%)" linearopts=(tickvaluepriority=true 
			 tickvaluesequence=(start=0 end=100 increment=10)) offsetmin=0.01 LABELATTRS=(Family= 'Courier New' Size=8pt))
  xaxisopts=(label='Time (day)' type=linear griddisplay=off display=(label tickvalues ticks) LABELATTRS=(Family= 'Courier New' Size=9)
            linearopts=(tickvaluepriority=true tickvaluelist=(&TimePointN.) 
			TICKVALUEFITPOLICY=NONE viewmin=0 viewmax=%sysfunc(ceil(&maxsur.)))
			offsetmin=0.02); 
         				    
  stepplot x=time y=SURVIVAL/group=grpn name="step";
  scatterplot x=time y=CNSR_SURVIVAL/group=grpn name="scatter";

  layout gridded/columns=3 border=false autoalign=auto columngutter=1;
    entry halign=left ""/TEXTATTRS=(size=8pt );	 entry halign=left "Treatment"/TEXTATTRS=(size=8pt Family= 'Courier New'); entry halign=left "Placebo"/TEXTATTRS=(size=8pt Family= 'Courier New');
    entry halign=left "Event/N"/TEXTATTRS=(size=8pt Family= 'Courier New');  entry halign=left "&eventN1."/TEXTATTRS=(size=8pt); entry halign=left "&eventN2."/TEXTATTRS=(size=8pt);  
    entry halign=left "Median AE-Free,(95% CI)"/TEXTATTRS=(size=8pt Family= 'Courier New'); entry halign=left "&medci1."/TEXTATTRS=(size=8pt);  entry halign=left "&medci2."/TEXTATTRS=(size=8pt); 
  endlayout;
  
 
  mergedlegend "step" "scatter"/title="Group " titleattrs=(size=9pt family="Courier New") location=outside 
                border=true order=rowmajor across=5 autoalign=(topright topright center) valueattrs=(family='Courier New' Size=9 ) ;
  endlayout;
						 
  layout overlay;
    entry halign=left "No. of Patients at Risk" /valign=top border=false TEXTATTRS=(size=8pt Family= 'Courier New');
  endlayout;
						 
  layout overlay / pad=(bottom=2% left=2% right=2%) border=false walldisplay=none 
		xaxisopts=(display=none linearopts=(tickvaluepriority=true tickvaluelist=(&TimePointN.))) x2axisopts=(display=none) y2axisopts=(display=none) 
		yaxisopts=(reverse=true type=discrete display=(tickvalues) tickvalueattrs=(size=7pt));
		
		scatterplot x=arxvalue y=grpn / group=grpn markercharacterattrs=(size=8pt) markercharacter=atrisk1 ;
  endlayout;

endlayout;
endgraph;
end;
run;


/*RTF format*/
ods escapechar="~";
ods path work.templat(update) sasuser.templat(read) sashelp.tmplmst(read);
options nocenter nonumber nobyline nodate nonumber;

ODS listing close;
ODS RTF file="C:\Users\feiy\OneDrive - The University of Colorado Denver\Desktop\ALEX\Real Data\Output\&out..rtf" style=custom  nogtitle nogfootnote;

ODS GRAPHICS ON/noborder outputfmt=png height=5in width=9in;

proc sgrender data=survial template=plot;
  format grpn grp.;
  run;

ODS RTF CLOSE;
ODS listing close;
ods graphics off;
%mend;

%kmrecov(in=az,out=az_km_ae);
%kmrecov(in=brii,out=brii_km_ae);
%kmrecov(in=lilly,out=lilly_km_ae);
%kmrecov(in=mp,out=mp_km_ae);
%kmrecov(in=pf,out=pf_km_ae);
%kmrecov(in=vir,out=vir_km_ae);
