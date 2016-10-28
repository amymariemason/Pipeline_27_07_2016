* Compare the three methods in phenotype predictions

set li 130

cap log close
log using analysis2.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

*********** create phenotype predictions by site
use anti_panel_all, clear

******************************************************
noi di _n(5) _dup(80) "=" _n " 1 phenotype predictions by site" _n _dup(80) "="

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

******************************************************
noi di _n(5) _dup(80) "=" _n " 2 phenotype predictions compared to gold standard" _n _dup(80) "="
use anti_panel_all, clear

* in table
noi tab gold valueall

* reduced table
drop if !inlist(gold, "R", "S")
noi tab gold valueall
gen agree = (valueall=="rrr" & gold=="R")|(valueall=="sss" & gold=="S")
summ agree
noi di r(sum) " out of " r(N) "results agree between all methods and gold standard"
noi di r(sum)/r(N)*100

* in graph

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




* sensitivity/ specificity

******************************************************
noi di _n(5) _dup(80) "=" _n " 3 find sensitivity and specifity of each method/antibiotic" _n _dup(80) "="

**************************
* TYPE WRITER
****************************

use anti_panel_all, clear

assert site!=""
assert _N==16548

noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "R", "S")
summ clear
noi di r(sum) " results are R/S"
drop if clear!=1


noi di "Typewriter"

contract  site  gold valuetype
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

* graph
use typewriter, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

#delimit ;
twoway rcap tw_usens tw_lsens num || scatter tw_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Typewriter vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_sens.tif", as(tif) replace


#delimit ;
twoway rcap tw_uspec tw_lspec num || scatter tw_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Typewriter vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_spec.tif", as(tif) replace


***************************************************
* GENEFINDER
***************************************************


use anti_panel_all, clear

assert site!=""
assert _N==16548

noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "R", "S")
summ clear
noi di r(sum) " results are R/S"
drop if clear!=1


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


* graph
use genefinder, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

#delimit ;
twoway rcap gf_usens gf_lsens num || scatter gf_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Genefinder vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\gf_sens.tif", as(tif) replace


#delimit ;
twoway rcap gf_uspec gf_lspec num || scatter gf_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Genefinder vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\gf_spec.tif", as(tif) replace



***************************************************
* MYKROBE
***************************************************


use anti_panel_all, clear

assert site!=""
assert _N==16548

noi di "restrict to sites with clear gold standard values only"

gen clear=inlist(gold, "R", "S")
summ clear
noi di r(sum) " results are R/S"
drop if clear!=1


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


* graph
use mykrobe, clear
sort site
gen num =_n 
labmask num, values(site)
local max= _N

#delimit ;
twoway rcap z_usens z_lsens num || scatter z_sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Mykrobe vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\z_sens.tif", as(tif) replace


#delimit ;
twoway rcap z_uspec z_lspec num || scatter z_spec num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Mykrobe vs. Lab") 
subtitle("Specificity in Phenotype prediction");
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\z_spec.tif", as(tif) replace

***************************
* COMBINED
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

#delimit ;
twoway rcap usens lsens newcount, ylabel(,format(%3.2f)) xlabel(4(5)59, valuelabel angle(90))
|| scatter  sens newcount if method=="gf", mcolor(orange)
||  scatter  sensitivity newcount if method =="tw", mcolor(red) 
|| scatter  sens newcount if method=="z", mcolor(blue)
legend(order(2 3 4) lab(2 "Genefinder") lab(3 "Typewriter") lab(4 "Mykrobe"))
subtitle("Sensitivity in Phenotype prediction") graphregion(fcolor(white));
#delimit cr

graph save Graph "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_sens.gph", replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_sens.tif", as(tif) replace

#delimit ;
twoway rcap uspec lspec newcount, ylabel(,format(%3.2f))
|| scatter  spec newcount if method=="gf", mcolor(orange)
|| scatter  spec newcount if method =="tw", mcolor(red)
|| scatter  spec newcount if method=="z", mcolor(blue)
xlabel(4(5)59, valuelabel angle(90))
legend(order(2 3 4) lab(2 "Genefinder") lab(3 "Typewriter") lab(4 "Mykrobe") )
subtitle("Specificity in Phenotype prediction") graphregion(fcolor(white));
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\all_spec.tif", as(tif) replace


**************************************************
* OVERALL
**************************************************




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
