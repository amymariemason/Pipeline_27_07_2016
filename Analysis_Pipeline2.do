************************************************
* ANALYSIS_PIPELINE2.DO 
************************************************

* Analyses the three methods on antibiotic predictions; creates summary tables and graphs

*Inputs: anti_panel_all  from create_predict_anti.do
* Outputs :  phenotype_disagreements.tif ( bar chart of method results by antibiotic)
*  phenotype_disagreements_withgold.tif ( bar charts of method results by antibiotic and lab results, restricted to results where lab result either R or S)
* all_sens.tif ( sensitivity of each method, by antibiotic) and all_spec.tif (specificity of each method, by antibiotic)
* graphs saved to Graph folder
* Written by: Amy Mason

***************************************************************************
* Compare the three methods in phenotype predictions for antibiotics
****************************************************************************
set li 130

cap log close
log using analysis2.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets


use anti_panel_all, clear

******************************************************
noi di _n(5) _dup(80) "=" _n " 1 phenotype predictions by site" _n _dup(80) "="
********************************************************
* create summary tables of the phenotype predictions
* note r/s = method predictions, R/S = lab results

noi tab site valuea, m
gen agree  = (valuea=="rrr"| valuea=="sss" )
summ agree
noi di r(sum) " out of " r(N) " phenotype predictions agree"
noi di r(sum)/_N*100

noi di "largest disagreements"
noi tab site if agree!=1, sort
noi tab site valuea if agree!=1

gen m_overcall = 1 if valuea=="srs"
summ m_overcall
noi di r(sum) " out of " r(N) "mykrobe overcall"
noi di r(sum)/_N*100

gen t_undercall = 1 if valuea=="rrs"
summ t_undercall
noi di r(sum) " out of " r(N) "typewriter undercall"
noi di r(sum)/_N*100


* create values for bar chart of all method results
noi di" make bar chart of which methods differ on which antibiotic"
contract  site  value*
rename valueall all
drop value*
reshape wide _freq, i(site) j(all) string
rename _freq* *


* bar graph of predictions (antibiotics)
graph bar (asis)  rrr rrs rsr rss srs sss,  over(site, label(angle(90))) stack title("Combinations of Results") subtitle("Results given as Genefinder Mykrobe Typewriter ") legend( label(1 "rrr") label( 2 "rrs" ) label (3 "rsr") label (4 "rss" ) label (5 "srs") label( 6 "sss" )) bar(1, color(gs8)) bar(2, color(cyan)) bar(3, color(red)) bar(4, color(blue)) bar(5, color(yellow)) bar(6, color(gs12))
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\phenotype_disagreements.tif", as(tif) width(2550) replace


*************************************************************
* Bar chart comparing with gold standard
*************************************************************
noi di _n(5) _dup(80) "=" _n " 2 phenotype predictions compared to gold standard" _n _dup(80) "="
use anti_panel_all, clear

* summary tables of results
noi tab gold valueall
gen agree = (valueall=="rrr")|(valueall=="sss")
summ agree
noi di r(sum) " out of " r(N) "results agree between all methods"
noi di r(sum)/r(N)*100


noi di "disagreements by site"
noi bysort gold: tab site valueall

noi di "discrenpancies"
sort site gold valueall sample site
noi list sample site valueall if agree!=1
drop agree


noi di "disagreements with gold standard"
drop if !inlist(gold, "R", "S")
noi tab gold valueall
gen agree = (valueall=="rrr" & gold=="R")|(valueall=="sss" & gold=="S")
summ agree
noi di r(sum) " out of " r(N) "results agree between all methods and gold standard"
noi di r(sum)/r(N)*100


noi di "disagreements by site"
noi bysort gold: tab site valueall

noi di "discrenpancies"
sort site gold valueall sample site
noi list sample site gold valueall if agree!=1


* creat bar graph values

contract  site  gold valueall

reshape wide _freq, i(site gold) j(valueall) string
rename _freq* *


sort gold site
 gen num =_n
 replace num=num+1 if gold=="S"
labmask num, values(site)
gen goldlabel = "Gold standard " + gold

* bar graph of predictions (antibiotics)
graph bar (asis) rrr rrs rsr rss srs sss if inlist(golds,"R", "S"),graphregion(color(white))    over(site, label(angle(90) labsize(tiny)) )  over(goldlabel, label(labsize(small)) ) stack  title("Prediction combinations") subtitle("Results given as Genefinder Mykrobe Typewriter") legend(rows(2) label(1 "rrr") label( 2 "rrs" ) label (3 "rsr") label (4 "rss" ) label (5 "srs") label( 6 "sss" )) bar(1, color(gs8)) bar(2, color(cyan)) bar(3, color(red)) bar(4, color(blue)) bar(5, color(yellow)) bar(6, color(gs12)) ylabel(0(500)1000 1400)
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\phenotype_disagreements_withgold.tif", as(tif) width(2550) replace



*********************************************************
* sensitivity/ specificity graphs for each method

******************************************************
noi di _n(5) _dup(80) "=" _n " 3 find sensitivity and specifity of each method/antibiotic" _n _dup(80) "="

**************************
* TYPE WRITER
****************************

use anti_panel_all, clear

assert site!=""
assert _N==16548

**
noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "R", "S")
summ clear
noi di r(sum) " results are R/S"
drop if clear!=1


noi di "Typewriter"

* contract data to get total values by site
contract  site  gold valuetype
rename value predict

assert predict!=""
* reshape the data to calculate the senstivity/ specificity

reshape wide _freq, i(site gold) j(predict) string
rename _freq* *

rename r combor
rename s combos
reshape wide combo*, i(site) j(gold) string
rename combo* *

* remember : resistance = postive result
rename sS TN
rename rR TP
rename sR FN 
*so predicted sensitive but actually resistant
rename rS FP 
* predicted resistant but actually sensitive

for any TN TP FN FP: replace X=0 if X==.

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

* save values for combining into multiple graphs

keep site  sens spec lsens usens lspec uspec
noi list  site  sens spec lsens usens lspec uspec
rename * tw_*
rename tw_site site
save typewriter, replace

* create labels
use typewriter, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

* make graph
#delimit ;
twoway rspike tw_usens tw_lsens num, lcolor(black) || scatter tw_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Typewriter vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_sens.tif", as(tif) replace



#delimit ;
twoway rspike tw_uspec tw_lspec num, lcolor(black) || scatter tw_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Typewriter vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr

* save graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_spec.tif", as(tif) replace


***************************************************
* GENEFINDER
***************************************************


use anti_panel_all, clear

assert site!=""
assert _N==16548

* restrict
noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "R", "S")
summ clear
noi di r(sum) " results are R/S"
drop if clear!=1

* contract to get totals for sites
noi di "Genefinder"
contract  site  gold valuegene
rename value predict
assert predict!=""

* reshape
reshape wide _freq, i(site gold) j(predict) string
rename _freq* *
rename r combor
rename s combos
reshape wide combo*, i(site) j(gold) string
rename combo* *

* remember : resistance = postive result
rename sS TN
rename rR TP
rename sR FN 
*so predicted sensitive but actually resistant
rename rS FP 
* predicted resistant but actually sensitive

for any TN TP FN FP: replace X=0 if X==.

* generate sens and spec
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
* make sens graph
#delimit ;
twoway rspike gf_usens gf_lsens num, lcolor(black) || scatter gf_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Genefinder vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
* save graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\gf_sens.tif", as(tif) replace

* make specificity graph
#delimit ;
twoway rspike gf_uspec gf_lspec num, lcolor(black) || scatter gf_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Genefinder vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr
* save graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\gf_spec.tif", as(tif) replace



***************************************************
* MYKROBE
***************************************************


use anti_panel_all, clear

assert site!=""
assert _N==16548

* restrict
noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "R", "S")
summ clear
noi di r(sum) " results are R/S"
drop if clear!=1

* contract
noi di "Mykrobe"
contract  site  gold valuezam
rename value predict
assert predict!=""

* reshape
reshape wide _freq, i(site gold) j(predict) string
rename _freq* *

rename r combor
rename s combos
reshape wide combo*, i(site) j(gold) string
rename combo* *

* remember : resistance = postive result
rename sS TN
rename rR TP
rename sR FN 
*so predicted sensitive but actually resistant
rename rS FP 
* predicted resistant but actually sensitive

for any TN TP FN FP: replace X=0 if X==.

* gen sens/ spec point values
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


* graph labels/ order
use mykrobe, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

* make sens graph
#delimit ;
twoway rspike z_usens z_lsens num, lcolor(black) || scatter z_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Mykrobe vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
* save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\z_sens.tif", as(tif) replace

* make spec graph
#delimit ;
twoway rspike z_uspec z_lspec num, lcolor(black) || scatter z_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Mykrobe vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr
* save
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\z_spec.tif", as(tif) replace

***************************
* COMBINE THE SENS/ SPEC GRAPHS
***************************

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
rename z_ meth_z

reshape long meth_, i(site new) j(method) string
rename meth_ value

reshape wide value, i(site method) j(new) string
rename value* *

order site method sens lsens usens spec lspec uspec


* sneakily creating gaps on graphs
expand 3 if method =="gf", gen(dups)
replace method ="blank" if dups==1
replace sens=. if method=="blank"
replace lsens=. if method=="blank"
replace usens=. if method=="blank"
replace uspec=. if method=="blank"
replace lspec=. if method=="blank"
replace spec=. if method=="blank"

* order the slabels
encode(method), gen(order)
sort site order
gen newcount = _n
local max =_N
replace site = "fusidic acid" if strpos(site, "fus")
labmask newcount, values(site)

* make sens graph
#delimit ;
twoway rspike usens lsens newcount, ylabel(,format(%3.2f)) xlabel(4(5)59, valuelabel angle(90)) lcolor(black)
|| scatter  sens newcount if method=="gf", mcolor(orange)
||  scatter  sensitivity newcount if method =="tw", mcolor(red) 
|| scatter  sens newcount if method=="z", mcolor(blue)
legend(order(2 3 4) lab(2 "Genefinder") lab(3 "Typewriter") lab(4 "Mykrobe") rows(1))
subtitle("Sensitivity in Phenotype prediction") graphregion(fcolor(white))
xtitle("Antibiotic");
#delimit cr

* save graph
graph save Graph "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_sens.gph", replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_sens.tif", as(tif) replace

* make spec graph
#delimit ;
twoway rspike uspec lspec newcount, ylabel(,format(%3.2f))  lcolor(black)
|| scatter  spec newcount if method=="gf", mcolor(orange)
|| scatter  spec newcount if method =="tw", mcolor(red)
|| scatter  spec newcount if method=="z", mcolor(blue)
xlabel(4(5)59, valuelabel angle(90))
legend(order(2 3 4) lab(2 "Genefinder") lab(3 "Typewriter") lab(4 "Mykrobe") rows(1))
subtitle("Specificity in Phenotype prediction") graphregion(fcolor(white))
xtitle("Antibiotic");
#delimit cr
* save graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_spec.tif", as(tif) replace


**************************************************
* OVERALL
**************************************************
* add the overall sensitivity and specificity to the log file



noi di "all methods"
tempfile temp
	use anti_panel_all, clear
	assert site!=""
	assert _N==16548
* restrict to clear goldstandard values
	gen clear=inlist(gold, "R", "S")
	summ clear
	noi di r(sum) " results are R/S"
	drop if clear!=1
	contract  gold valuegene
	rename value predict
	assert predict!=""
* reshape
	reshape wide _freq, i( gold) j(predict) string
	rename _freq* *
	rename r combor
	rename s combos
	gen method = "valuegene"
	reshape wide combo*, i(method) j(gold) string
	rename combo* *

	* remember : resistance = postive result
	rename sS TN
	rename rR TP
	rename sR FN 
	*so predicted sensitive but actually resistant
	rename rS FP 
	* predicted resistant but actually sensitive
save temp, replace
	
foreach k in valuetype valuez{
	use anti_panel_all, clear
	assert site!=""
	assert _N==16548
* restrict to clear goldstandard values
	gen clear=inlist(gold, "R", "S")
	summ clear
	noi di r(sum) " results are R/S"
	drop if clear!=1
	contract  gold `k'
	rename value predict
	assert predict!=""
* reshape
	reshape wide _freq, i( gold) j(predict) string
	rename _freq* *
	rename r combor
	rename s combos
	gen method = "`k'"
	reshape wide combo*, i(method) j(gold) string
	rename combo* *

	* remember : resistance = postive result
	rename sS TN
	rename rR TP
	rename sR FN 
	*so predicted sensitive but actually resistant
	rename rS FP 
	* predicted resistant but actually sensitive
	append using temp
	save temp, replace
}

	
* point estimates
	for any TN TP FN FP: replace X=0 if X==.

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


exit

* after this point was written to compare which sites where causing the discrepencies when mykrobe was acting up
* not relevant to paper

************************
*antibiotic specific: rifampicin, fusicidic acid, pencilin
***********************
* rpob
noi di "results for rifamcipin (site: rpob)"
use anti_panel_all, clear
keep if strpos(site, "rif")
noi di "Recall order is  Typewriter Mykrobe Genefinder (lower case) Gold (uppercase)"
noi tab gold valuea

*fusa, fusb, fusc
noi di "results for fusidic acid: fusb, fusc, fusa"
use anti_panel_all, clear
keep if strpos(site, "fus")
noi di "Recall order is  Typewriter Mykrobe Genefinder (lower case) Gold (uppercase)"
noi tab gold valuea
* split into site by site results
tempfile fustemp
keep sample gold valuez
save fustemp, replace

use pipeline_clean_all_values_wide, clear
keep if inlist(site, "fusa", "fusb", "fusc")
noi tab site ValueA
keep sample site ValueA
reshape wide ValueA, i(sample) j(site) string
merge 1:1 sample using fustemp, update
assert _merge==3
drop _merge

gen ABC= "A:" + ValueAllfusa + ",B:"+  ValueAllfusb + ",C:" + ValueAllfusc
noi tab ABC gold
keep if gold=="R"
noi di "for resistant samples"
noi tab ABC valuez


* meca, mecc, blaz
noi di "results for penicillin: meca, mecc, blaz"
use anti_panel_all, clear
keep if strpos(site, "pen")
noi di "Recall order is  Typewriter Mykrobe Genefinder (lower case) Gold (uppercase)"
noi tab gold valuea
* split into site by site results
tempfile pentemp
keep sample gold valuez
save pentemp, replace

use pipeline_clean_all_values_wide, clear
keep if inlist(site, "meca", "mecc", "blaz")
noi tab site ValueA
keep sample site ValueA
reshape wide ValueA, i(sample) j(site) string
merge 1:1 sample using pentemp, update
assert _merge==3
drop _merge

gen ABC= "A:" + ValueAllmeca + ",C:"+  ValueAllmecc + ",Z:" + ValueAllblaz
noi tab ABC gold
keep if gold=="R"
noi di "for resistant samples"
noi tab ABC valuez


*cipro
noi di "results for ciprofloxin: fusb, fusc, fusa"
use anti_panel_all, clear
keep if strpos(site, "cipro")
noi di "Recall order is  Typewriter Mykrobe Genefinder (lower case) Gold (uppercase)"
noi tab gold valuea
* split into site by site results
tempfile ciprotemp
keep sample gold valuez
save ciprotemp, replace

use pipeline_clean_all_values_wide, clear
keep if inlist(site, "grla/gyra")
noi tab site ValueA
keep sample site ValueA
merge 1:1 sample using ciprotemp, update
assert _merge==3
drop _merge

noi tab ValueAll gold
keep if gold=="R"
noi di "for resistant samples"
noi tab ValueAll valuez



* meca, mecc, blaz
noi di "results for erthytomycin: msta, mphc, erma, ermb, ermc, ermt, ermy"
use anti_panel_all, clear
keep if strpos(site, "ery")
noi di "Recall order is  Typewriter Mykrobe Genefinder (lower case) Gold (uppercase)"
noi tab gold valuea
* split into site by site results
tempfile ertemp
keep sample gold valuez
save ertemp, replace

use pipeline_clean_all_values_wide, clear
keep if inlist(site, "erma", "ermb", "ermc", "ermt", "ermy", "msra", "mphc" )
noi tab site ValueA
keep sample site ValueA
reshape wide ValueA, i(sample) j(site) string
merge 1:1 sample using ertemp, update
assert _merge==3
drop _merge

gen ABC= "a:" + ValueAllerma + ",B:"+  ValueAllermb + ",C:" + ValueAllermc + ",T:" + ValueAllermt + ",Y:" + ValueAllermy + ",M1:" + ValueAllmsra + ",M2:" + ValueAllmphc
noi tab ABC gold
noi di "drop those with complete agreement on all sites and with gold"
gen P_present = strpos(ABC, "P")
drop if P_present==0 & gold=="S"
gen P_agree = strpos(ABC, "PPP")
drop if P_agree!=0 & gold=="R"
sort ABC sample

noi di "resistant gold"
noi list sample ABC if gold=="R"

noi di "sensitive gold"
noi list sample ABC if gold=="S"