************************************************
* ANALYSIS_PIPELINE3.DO 
************************************************

* Analyses the three methods on virulence factor predictions; creates summary tables and graphs

*Inputs: viru_panel_all  from create_predict_viru.do
* Outputs : virulence_disagreements.wmf ( bar chart of method results by antibiotic)
*  virulence_disagreements_withgold.wmf ( bar charts of method results by antibiotic and lab results, restricted to results where lab result either R or S)
* all_viru_sens.wmf ( sensitivity of each method, by virulence factor) and all_viru_spec.wmf (specificity of each method, by virulence factor)
* graphs saved to Graph folder
* Written by: Amy Mason

set li 130

cap log close
log using analysis3.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

*********************************************
use viru_panel_all, clear

******************************************************
noi di _n(5) _dup(80) "=" _n " 1 virulence predictions by site" _n _dup(80) "="
* create summary tables

noi tab site valuea, m
gen agree  = (valuea=="ppp"| valuea=="aaa" )
summ agree
noi di r(sum) " out of " r(N) " virulence predictions agree"
noi di r(sum)/_N*100


noi di "disagreements by site"
noi bysort gold: tab site valueall

noi di "discrenpancies"
sort site gold valueall sample site
noi list sample site valueall if agree!=1



noi di "largest disagreements"
noi tab site if agree!=1, sort
noi tab site valuea if agree!=1


gen t_undercall = (valuea=="ppa")
summ t_undercall
noi di r(sum) " out of " r(N) "typewriter undercall"
noi di r(sum)/_N*100



* contract to make first bar chart
noi di" make bar chart of which methods differ on which virulence factor"
preserve
contract  site  value*
rename valueall all
drop value*

reshape wide _freq, i(site) j(all) string
rename _freq* *


* bar graph of predictions (antibiotics)
#delimit ;
graph bar (asis) ppp app ppa apa aaa,  over(site, label(angle(90))) stack title("Combinations of Results") subtitle("Results given as Genefinder Mykrobe Typewriter ") 
legend( label(1 "ppp") label( 2 "app" ) label (3 "ppa") label (4 "apa" ) label (5 "aaa") ) 
bar(1, color(gs8)) bar(2, color(green)) bar(3, color(cyan)) bar(4, color(yellow)) bar(5, color(gs12));
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\virulence_disagreements.wmf", as(wmf) replace

restore

drop agree

*************************************************************
*  create bar chart comparing virulence to lab results (restricted to those we have results for)

******************************************************
noi di _n(5) _dup(80) "=" _n " 2 phenotype predictions compared to gold standard" _n _dup(80) "="


* in table
noi di "overall results"
noi tab gold valueall

noi di "results by site"
noi tab site valueall


* reduced table
noi di "comparision to clear gold value results"
drop if !inlist(gold, "A", "P")
noi tab gold valueall
gen agree = (valueall=="ppp" & gold=="P")|(valueall=="aaa" & gold=="A")
summ agree
noi di r(sum) " out of " r(N) "results agree between all methods and gold standard"
noi di r(sum)/r(N)*100

noi di "disagreements by site"
noi bysort gold: tab site valueall

* contract to make graph values
preserve
contract  site  gold valueall

reshape wide _freq, i(site gold) j(valueall) string
rename _freq* *


gen graphx = site + " " + gold
drop if !inlist(gold, "A", "P")
sort gold site
 gen num =_n
 replace num=num+1 if gold=="P"
labmask num, values(site)
gen goldlabel = "Gold standard " + gold

* bar graph of predictions  compared to lab results
#delimit ;
graph bar (asis) ppp  ppa apa aaa if inlist(golds,"A", "P"),graphregion(color(white))    
over(site, label(angle(90) labsize(tiny)) )  over(goldlabel, label(labsize(small)) ) stack  
ytitle("Number of isolates")  
legend( label(1 "ppp")  label (2 "ppa") label (3 "apa" ) label (4 "aaa") subtitle ("Results given as Genefinder Mykrobe Typewriter") ) 
bar(1, color(gs8))  bar(2, color(cyan)) bar(3, color(yellow)) bar(4, color(gs12))
ylabel(0(50)250);
#delimit cr
*save graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\virulence_disagreements_withgold.wmf", as(wmf) replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\virulence_disagreements_withgold.tif", width(2280) replace

restore
*********************************************************************
* sensitivity/ specificity for each method/ virulance factor

******************************************************
noi di _n(5) _dup(80) "=" _n " 3 find sensitivity and specifity of each method/virulence" _n _dup(80) "="

**************************
* TYPE WRITER
****************************

use viru_panel_all, clear

assert site!=""
assert _N==26201

* drop samples without lab results
noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "P", "A")
summ clear
noi di r(sum) " results are P/A"
drop if clear!=1

* contract to get result totals
noi di "Typewriter"

contract  site  gold valuetype
rename value predict

assert predict!=""

* reshape

reshape wide _freq, i(site gold) j(predict) string
rename _freq* *

rename p combop
rename a comboa
reshape wide combo*, i(site) j(gold) string
rename combo* *

* remember : resistance = postive result
rename aA TN
rename pP TP
rename aP FN 
*so predicted sensitive but actually resistant
rename pA FP 
* predicted resistant but actually sensitive

for any TN TP FN FP: replace X=0 if X==.

* point estimates for sensitivity and specificity 
gen sensitivity = TP/(TP+FN)
gen specificity = TN/(TN+FP)

gen predictPOS = TP+FP
gen predictNEG = TN +FN

gen trueNeg= TN +FP
gen truePos= TP +FN

gen MajorErrorRate = FN/truePos
gen VeryMajorErrorRate = FP/trueNeg
gen Agreement = (TP+TN)/(trueNeg+truePos)

* get confidence intervals
for any lsens usens lspec uspec lme ume lvme uvme:gen X=.
local max=_N
forvalues i=1(1)`max'{
if truePos[`i']>0{
cii prop  truePos[`i'] TP[`i']
replace lsens = r(lb) if _n==`i'
replace usens = r(ub) if _n==`i'
cii prop truePos[`i'] FN[`i']
replace lvme = r(lb) if _n==`i'
replace uvme = r(ub) if _n==`i'
}
else{
noi di site[`i'] " has no positive samples"
}
if trueNeg[`i']>0{
cii prop trueNeg[`i'] TN[`i']
replace lspec = r(lb) if _n==`i'
replace uspec = r(ub) if _n==`i'
cii prop trueNeg[`i'] FP[`i']
replace lme = r(lb) if _n==`i'
replace ume = r(ub) if _n==`i'
}
else{
noi di site[`i'] " has no negative samples"
}
}

* save values for combining into graphs


keep site  sens spec lsens usens lspec uspec
noi list  site  sens spec lsens usens lspec uspec
rename * tw_*
rename tw_site site
save typewriter, replace

* graph ordering/labels
use typewriter, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

* make sens graph
#delimit ;
twoway rspike tw_usens tw_lsens num, lcolor(black) || scatter tw_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Typewriter vs. Lab") 
subtitle("Sensitivity in Virulence prediction");
#delimit cr
*save graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_viru_sens.wmf", as(wmf) replace

* make spec graph
#delimit ;
twoway rspike tw_uspec tw_lspec num, lcolor(black) || scatter tw_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Typewriter vs. Lab") 
subtitle("Specificity in Virulence prediction");
#delimit cr
* save graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_viru_spec.wmf", as(wmf) replace


***************************************************
* GENEFINDER
***************************************************
use viru_panel_all, clear

assert site!=""
assert _N==26201
* restrict
noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "A", "P")
summ clear
noi di r(sum) " results are P/A"
drop if clear!=1

* contract
noi di "Genefinder"

contract  site  gold valuegene
rename value predict

assert predict!=""
* reshape

reshape wide _freq, i(site gold) j(predict) string
rename _freq* *

rename p combop
rename a comboa
reshape wide combo*, i(site) j(gold) string
rename combo* *

* remember : resistance = postive result
rename aA TN
rename pP TP
rename aP FN 
*so predicted sensitive but actually resistant
rename pA FP 
* predicted resistant but actually sensitive

for any TN TP FN FP: replace X=0 if X==.
* point estimates
gen sensitivity = TP/(TP+FN)
gen specificity = TN/(TN+FP)

gen predictPOS = TP+FP
gen predictNEG = TN +FN

gen trueNeg= TN +FP
gen truePos= TP +FN

gen MajorErrorRate = FN/truePos
gen VeryMajorErrorRate = FP/trueNeg
gen Agreement = (TP+TN)/(trueNeg+truePos)

* get confidence intervals
for any lsens usens lspec uspec lme ume lvme uvme:gen X=.
local max=_N
forvalues i=1(1)`max'{
if truePos[`i']>0{
cii prop  truePos[`i'] TP[`i']
replace lsens = r(lb) if _n==`i'
replace usens = r(ub) if _n==`i'
cii prop truePos[`i'] FN[`i']
replace lvme = r(lb) if _n==`i'
replace uvme = r(ub) if _n==`i'
}
else{
noi di site[`i'] " has no positive samples"
}
if trueNeg[`i']>0{
cii prop trueNeg[`i'] TN[`i']
replace lspec = r(lb) if _n==`i'
replace uspec = r(ub) if _n==`i'
cii prop trueNeg[`i'] FP[`i']
replace lme = r(lb) if _n==`i'
replace ume = r(ub) if _n==`i'
}
else{
noi di site[`i'] " has no negative samples"
}
}

* save values for combining into graphs


keep site  sens spec lsens usens lspec uspec
noi list  site  sens spec lsens usens lspec uspec
rename * gf_*
rename gf_site site
save genefinder, replace


* graph labels
use genefinder, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

* sensitivity graph
#delimit ;
twoway rspike gf_usens gf_lsens num, lcolor(black) || scatter gf_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Genefinder vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
* save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\gf_viru_sens.wmf", as(wmf) replace

* specificity graph
#delimit ;
twoway rspike gf_uspec gf_lspec num, lcolor(black) || scatter gf_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Genefinder vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr
* save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\gf_viru_spec.wmf", as(wmf) replace



***************************************************
* MYKROBE
***************************************************
use viru_panel_all, clear

assert site!=""
assert _N==26201
* restrict
noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "A", "P")
summ clear
noi di r(sum) " results are P/A"
drop if clear!=1

* contract
noi di "Mykrobe"

contract  site  gold valuezam
rename value predict

assert predict!=""
* reshape

reshape wide _freq, i(site gold) j(predict) string
rename _freq* *

rename p combop
rename a comboa
reshape wide combo*, i(site) j(gold) string
rename combo* *

* remember : resistance = postive result
rename aA TN
rename pP TP
rename aP FN 
*so predicted sensitive but actually resistant
rename pA FP 
* predicted resistant but actually sensitive

for any TN TP FN FP: replace X=0 if X==.
* point estimates
gen sensitivity = TP/(TP+FN)
gen specificity = TN/(TN+FP)

gen predictPOS = TP+FP
gen predictNEG = TN +FN

gen trueNeg= TN +FP
gen truePos= TP +FN

gen MajorErrorRate = FN/truePos
gen VeryMajorErrorRate = FP/trueNeg
gen Agreement = (TP+TN)/(trueNeg+truePos)

* get confidence intervals
for any lsens usens lspec uspec lme ume lvme uvme:gen X=.
local max=_N
forvalues i=1(1)`max'{
if truePos[`i']>0{
cii prop  truePos[`i'] TP[`i']
replace lsens = r(lb) if _n==`i'
replace usens = r(ub) if _n==`i'
cii prop truePos[`i'] FN[`i']
replace lvme = r(lb) if _n==`i'
replace uvme = r(ub) if _n==`i'
}
else{
noi di site[`i'] " has no positive samples"
}
if trueNeg[`i']>0{
cii prop trueNeg[`i'] TN[`i']
replace lspec = r(lb) if _n==`i'
replace uspec = r(ub) if _n==`i'
cii prop trueNeg[`i'] FP[`i']
replace lme = r(lb) if _n==`i'
replace ume = r(ub) if _n==`i'
}
else{
noi di site[`i'] " has no negative samples"
}
}

* save values for combining into graphs


keep site  sens spec lsens usens lspec uspec
noi list  site  sens spec lsens usens lspec uspec
rename * z_*
rename z_site site
save mykrobe, replace


* graph labels
use mykrobe, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

* sens graph
#delimit ;
twoway rspike z_usens z_lsens num, lcolor(black) || scatter z_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Mykrobe vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
*save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\z_viru_sens.wmf", as(wmf) replace

*spec graph
#delimit ;
twoway rspike z_uspec z_lspec num, lcolor(black) || scatter z_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Mykrobe vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr
*save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\z_viru_spec.wmf", as(wmf) replace


********************************************
* COMBINE THE SENSITIVITY AND SPECIFICITY GRAPHS
*******************************************

* merge the three sets
use typewriter, clear
 merge 1:1 site using genefinder, update
 assert _merge==3
 drop _merge
 merge 1:1 site using mykrobe, update 
 assert _merge==3
 drop _merge
 
 * reshape to enable graph
reshape long tw_ gf_ z_, i(site) j(new) string
rename tw_ meth_tw
rename gf_ meth_gf
rename z_ meth_m

reshape long meth_, i(site new) j(method) string
rename meth_ value

reshape wide value, i(site method) j(new) string
rename value* *

order site method sens lsens usens spec lspec uspec


* sneaky creating gaps on graphs
expand 3 if method =="gf", gen(dups)
replace method ="blank" if dups==1
replace sens=. if method=="blank"
replace lsens=. if method=="blank"
replace usens=. if method=="blank"
replace uspec=. if method=="blank"
replace lspec=. if method=="blank"
replace spec=. if method=="blank"

* order bits
encode(method), gen(order)
sort site order
gen newcount = _n
local max =_N
replace site = "fusidic acid" if strpos(site, "fus")
labmask newcount, values(site)

* sensitivity grtaph
#delimit ;
twoway rspike usens lsens newcount, ylabel(,format(%3.2f)) xlabel(4(5)94, valuelabel angle(90)) lcolor(black)
|| scatter  sens newcount if method=="gf", mcolor(orange) msize(small)
|| scatter  sens newcount if method=="m", mcolor(blue) msize(small)
||  scatter  sensitivity newcount if method =="tw", mcolor(red)  msize(small)
legend(order(2 3 4) lab(2 "Genefinder") lab(3 "Mykrobe") lab(4 "Typewriter")rows(1))
graphregion(fcolor(white))
ytitle("Sensitivity")
xtitle("Virulence Factor")
yscale(range(0.2 1))
ylabel(0.2 (0.2) 1);
#delimit cr
*save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_viru_sens.wmf", as(wmf) replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_viru_sens.tif", width(2880) replace

* specificity graph
#delimit ;
twoway rspike uspec lspec newcount, ylabel(,format(%3.2f)) lcolor(black)
|| scatter  spec newcount if method=="gf", mcolor(orange) msize(small)
|| scatter  spec newcount if method=="m", mcolor(blue) msize(small)
|| scatter  spec newcount if method =="tw", mcolor(red) msize(small)

xlabel(4(5)94, valuelabel angle(90))
legend(order(2 3 4) lab(2 "Genefinder") lab(3 "Mykrobe") lab(4 "Typewriter") rows(1))
ytitle("Specificity") graphregion(fcolor(white))
xtitle("Virulence Factor")
yscale(range(0.2 1))
ylabel(0.2 (0.2) 1)
;
#delimit cr
*save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_viru_spec.wmf", as(wmf) replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_viru_spec.tif", width(2880) replace

**************************************************
* OVERALL SPEC/SENS SUMMARY
**************************************************
* add overall sens/spec values to the log file


noi di "all methods"
tempfile temp
	use viru_panel_all, clear
	assert site!=""
	assert _N==26201
* restrict to clear goldstandard values
	gen clear=inlist(gold, "A", "P")
	summ clear
	noi di r(sum) " results are A/P"
	drop if clear!=1
	contract  gold valuegene
	rename value predict
	assert predict!=""
* reshape
	reshape wide _freq, i( gold) j(predict) string
	rename _freq* *
	rename p combop
	rename a comboa
	gen method = "valuegene"
	reshape wide combo*, i(method) j(gold) string
	rename combo* *

* remember : resistance = postive result
rename aA TN
rename pP TP
rename aP FN 
*so predicted sensitive but actually resistant
rename pA FP 
* predicted resistant but actually sensitive
save temp, replace
	
foreach k in valuetype valuez{
	use viru_panel_all, clear
	assert site!=""
	assert _N==26201
* restrict to clear goldstandard values
	gen clear=inlist(gold, "A", "P")
	summ clear
	noi di r(sum) " results are A/P"
	drop if clear!=1
	contract  gold `k'
	rename value predict
	assert predict!=""
* reshape
	reshape wide _freq, i( gold) j(predict) string
	rename _freq* *
	rename p combop
	rename a comboa
	gen method = "`k'"
	reshape wide combo*, i(method) j(gold) string
	rename combo* *

	* remember : resistance = postive result
	rename aA TN
	rename pP TP
	rename aP FN 
	*so predicted sensitive but actually resistant
	rename pA FP 
	* predicted resistant but actually sensitive
	append using temp
	save temp, replace
}

	for any TN TP FN FP: replace X=0 if X==.
* point estimates
	gen sensitivity = TP/(TP+FN)
	gen specificity = TN/(TN+FP)

	gen predictPOS = TP+FP
	gen predictNEG = TN +FN

	gen trueNeg= TN +FP
	gen truePos= TP +FN

	gen MajorErrorRate = FN/truePos
	gen VeryMajorErrorRate = FP/trueNeg
	gen Agreement = (TP+TN)/(trueNeg+truePos)

	* get confidence intervals
	for any lsens usens lspec uspec lme ume lvme uvme:gen X=.
	local max=_N
	forvalues i=1(1)`max'{
	if truePos[`i']>0{
	cii prop  truePos[`i'] TP[`i']
	replace lsens = r(lb) if _n==`i'
	replace usens = r(ub) if _n==`i'
	cii prop truePos[`i'] FN[`i']
	replace lvme = r(lb) if _n==`i'
	replace uvme = r(ub) if _n==`i'
	}
	else{
	noi di site[`i'] " has no positive samples"
	}
	if trueNeg[`i']>0{
	cii prop trueNeg[`i'] TN[`i']
	replace lspec = r(lb) if _n==`i'
	replace uspec = r(ub) if _n==`i'
	cii prop trueNeg[`i'] FP[`i']
	replace lme = r(lb) if _n==`i'
	replace ume = r(ub) if _n==`i'
	}
	else{
	noi di site[`i'] " has no negative samples"
	}
	}

* save values for combining for table
noi list method sens lsens usens spec lspec uspec
