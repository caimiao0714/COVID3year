/* This macro conducts generalized estimating equations using a log link and a Poisson distribution.
The events and offset term were both weighted,
95% CIs were generated from 2.5th and 97.5th percentiles of parametric bootstrapping 
The by-year estimates were conducted at day 359, 719, and 1079.
*/

%MACRO create_i(i);
%IF &i < 360 %THEN %LET nmax = 2;
%ELSE %IF 360 <= &i AND &i < 720 %THEN %LET nmax = 3;
%ELSE %IF 720 <= &i %THEN %LET nmax = 4;

PROC SQL;
CREATE TABLE dt&i AS
SELECT R.ScrSSN, R2.group AS group1, &i AS date, R2.lntime, 
	R.cvd_&i AS cvd, R.coag_&i AS coag, R.gi_&i AS gi, R.general_&i AS general, 
	R.kidney_&i AS kidney, R.mental_&i AS mental, R.muscol_&i AS muscol, R.neuro_&i AS neuro, 
	R.pulm_&i AS pulm, R.metab_&i AS metab, R.pasc15_&i AS pasc
	FROM sum20.sum20_&i AS R
LEFT JOIN ctem.DALY_mar_lntime_tvhc AS R2
	ON R.ScrSSN = R2.ScrSSN;
QUIT;
%MEND;




%MACRO GEE_rate(og_outcome, dateleft, dateright);
ods output GEEEmpPEst=pe ;
ods output GEENcov=cov;
proc genmod data=GEE_FULL ;
class ScrSSN group(ref="1 Control");
model &og_outcome = group period group*period / dist=poisson offset=lntime;
repeated subject=ScrSSN / mcovb;
run;


data pe1;
set  pe;
length _name_ $128.;
retain intercept  cov1 cov2 cov3 cov4 cov5;
_type_='PARMS';
_name_='time';
if _n_=1 then intercept=estimate;
if _n_=2 then cov1=estimate;
if _n_=3 then cov2=estimate;
if _n_=5 then cov3=estimate;
if _n_=6 then cov4=estimate;
if _n_=7 then cov5=estimate;
if _n_=7 then output;
keep intercept  cov1 cov2 cov3 cov4 cov5 _type_ _name_;
run;

%PUT "---------- Step 3 finish: PASC 2.0 GEE modeling; outcome: (&og_outcome)/&list ----------";
data cov1;
set cov;
length _name_ $128.;
_type_='COV';
if _n_=1 then _NAME_='intercept';
if _n_=2 then _NAME_='cov1';
if _n_=3 then _NAME_='cov2';
if _n_=4 then _NAME_='cov3';
if _n_=5 then _NAME_='cov4';
if _n_=6 then _NAME_='cov5';
if _name_^="";
run;


proc iml;
use cov1;
read all var _all_ into cov2;
use pe1;
read all var _all_ into pe2;
Numsamples=10000;
call randseed(12);
X=randnormal(numsamples,pe2,cov2);
create zzz from x;
append from x;
close zzz;
quit;

data zzz;
set  zzz;
rnum=_n_;
z=1;
rename col1=bint col2=bcov1 col3=bcov2 col4=bcov3 col5=bcov4 col6=bcov5;
run;


data hr;
set  zzz;

r00=exp(bint)*1050;
r01=exp(bint + bcov3)*1050;
r10=exp(bint + bcov1)*1050;
r11=exp(bint + bcov1 + bcov3 + bcov4)*1050;
r20=exp(bint + bcov2)*1050;
r21=exp(bint + bcov2 + bcov3 + bcov5)*1050;

diff0100 = r01 - r00;
diff1110 = r11 - r10;
diff2120 = r21 - r20;

irr1 = diff1110/diff0100;
irr2 = diff2120/diff0100;

cirr10 = r10/r00;
cirr20 = r20/r00;
cirr11 = r11/r01;
cirr21 = r21/r01;

burd1 = diff1110 - diff0100;
burd2 = diff2120 - diff0100;

diff1000 = r10 - r00;
diff2000 = r20 - r00;
diff2010 = r20 - r10;

diff1101 = r11 - r01;
diff2101 = r21 - r01;
diff2111 = r21 - r11;


num=_n_;
LENGTH outcome $16.;
outcome="&og_outcome";
run;



proc univariate data=hr;
var r00 r01 r10 r11 r20 r21 diff0100 diff1110 diff2120 irr1 irr2 
	cirr10 cirr20 cirr11 cirr21 
	burd1 burd2 diff1000 diff2000 diff2010 diff1101 diff2101 diff2111;
output out=hrci pctlpts=2.5 50 97.5 
pctlpre=r00_ r01_ r10_ r11_ r20_ r21_ diff0100_ diff1110_ diff2120_ irr1_ irr2_
	cirr10_ cirr20_ cirr11_ cirr21_ 
	burd1_ burd2_ diff1000_ diff2000_ diff2010_ diff1101_ diff2101_ diff2111_;
run;

data PASC20_&og_outcome;
set  hrci;
LENGTH outcome $16.;
DALYwt = "PASC2.0";
outcome = "&og_outcome";
date0 = &dateleft;
date1 = &dateright;

irr1 = catx("",round(irr1_50,0.01),"(",round(irr1_2_5,0.01),",",round(irr1_97_5,0.01),")");
irr2 = catx("",round(irr2_50,0.01),"(",round(irr2_2_5,0.01),",",round(irr2_97_5,0.01),")");

cirr10 = catx("",round(cirr10_50,0.01),"(",round(cirr10_2_5,0.01),",",round(cirr10_97_5,0.01),")");
cirr20 = catx("",round(cirr20_50,0.01),"(",round(cirr20_2_5,0.01),",",round(cirr20_97_5,0.01),")"); 
cirr11 = catx("",round(cirr11_50,0.01),"(",round(cirr11_2_5,0.01),",",round(cirr11_97_5,0.01),")"); 
cirr21 = catx("",round(cirr21_50,0.01),"(",round(cirr21_2_5,0.01),",",round(cirr21_97_5,0.01),")"); 

burd1 = catx("",round(burd1_50,0.01),"(",round(burd1_2_5,0.01),",",round(burd1_97_5,0.01),")");
burd2 = catx("",round(burd2_50,0.01),"(",round(burd2_2_5,0.01),",",round(burd2_97_5,0.01),")");

r00 = catx("",round(r00_50,0.01),"(",round(r00_2_5,0.01),",",round(r00_97_5,0.01),")");
r01 = catx("",round(r01_50,0.01),"(",round(r01_2_5,0.01),",",round(r01_97_5,0.01),")");
r10 = catx("",round(r10_50,0.01),"(",round(r10_2_5,0.01),",",round(r10_97_5,0.01),")");
r11 = catx("",round(r11_50,0.01),"(",round(r11_2_5,0.01),",",round(r11_97_5,0.01),")");
r20 = catx("",round(r20_50,0.01),"(",round(r20_2_5,0.01),",",round(r20_97_5,0.01),")");
r21 = catx("",round(r21_50,0.01),"(",round(r21_2_5,0.01),",",round(r21_97_5,0.01),")");

diff1000 = catx("",round(diff1000_50,0.01),"(",round(diff1000_2_5,0.01),",",round(diff1000_97_5,0.01),")");
diff2000 = catx("",round(diff2000_50,0.01),"(",round(diff2000_2_5,0.01),",",round(diff2000_97_5,0.01),")");
diff2010 = catx("",round(diff2010_50,0.01),"(",round(diff2010_2_5,0.01),",",round(diff2010_97_5,0.01),")");

diff1101 = catx("",round(diff1101_50,0.01),"(",round(diff1101_2_5,0.01),",",round(diff1101_97_5,0.01),")");
diff2101 = catx("",round(diff2101_50,0.01),"(",round(diff2101_2_5,0.01),",",round(diff2101_97_5,0.01),")");
diff2111 = catx("",round(diff2111_50,0.01),"(",round(diff2111_2_5,0.01),",",round(diff2111_97_5,0.01),")");

keep outcome DALYwt date0 date1 irr1 irr2 burd1 burd2 
	 cirr10 cirr20 cirr11 cirr21 r00 r01 r10 r11 r20 r21 
	 diff1000 diff2000 diff2010 diff1101 diff2101 diff2111;
run;
%MEND;




%let list=
pasc
kidney
cvd
coag
gi
general
mental
muscol
neuro
pulm
metab
;

%MACRO list_PASC20;
%DO j = 1 %TO 11;
		%let list_i=%scan(&list, &j);
		 PASC20_&list_i
%END;
%MEND;


%MACRO GEE_bootstrap(left, right);
%PUT "---------- Step 1A. create data &left ----------";
%create_i(&left);
%PUT "---------- Step 1B. create data &right ----------";
%create_i(&right);

%PUT "---------- Step 1C. Combine left and right data ----------";
DATA GEE_FULL;
SET  dt&left dt&right;
	LENGTH group $32.;
	IF date = &left THEN period = 0;
	ELSE period = 1;
	IF group1 = 'Control' THEN group = "1 Control";
	ELSE IF group1 = 'Covid nonhosp' THEN group = "2 Covid nonhosp";
	ELSE IF group1 = 'Covid hosp' THEN group = "3 Covid hosp";
RUN;

%PUT "========== Step 2. PASC 2.0 Modeling: start=&left, end: &right ==========";
%do q = 1 %to 11;
%let manifest=%scan(&&list,&q);
%PUT "---------- Step 2. &q/11. PASC 2.0 GEE modeling; outcome: &&manifest; start=&left, end: &right ----------";
%GEE_rate(&&manifest, &left, &right);

%end;

%PUT "========== Step 3. PASC 2.0 Combining organ-systems: start=&left, end: &right ==========";
DATA DALYpath.PASC20_&left._&right.;
SET  %list_PASC20;
RUN;

%MEND;
