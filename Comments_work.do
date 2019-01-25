************************************************
* COMMENTS_WORK.DO 
************************************************

*Additional tables etc as requested in feedback

*Inputs: 
* Outputs :  

* Written by: Amy Mason

***************************************************************************
* 
****************************************************************************
set li 130

cap log close
log using comments.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Datasets"


use anti_panel_all, clear

******************************************************
noi di _n(5) _dup(80) "=" _n " 1 prediction panel" _n _dup(80) "="
********************************************************

*antibiotics
use anti_panel_all, clear

gen all_value = valueg + valuez + valuet+ ": "+gold
replace all_value=valueg + ": "+gold if valueg==valuez & valuez==valuet
keep sample site set all
reshape wide all_value, i(sample set) j(site) string
rename all_value* *

save predictions, replace

* virulence

use viru_panel_all, clear

gen all_value = valueg + valuez + valuet+ ": "+gold
replace all_value=valueg+ ": "+gold if valueg==valuez & valuez==valuet
keep sample site set all
reshape wide all_value, i(sample set) j(site) string
rename all_value* *

merge 1:1 sample using predictions, update
assert _merge==3
drop _merge

* export
export excel using "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Comments\Predictions.xls", sheet("Predictions") firstrow(variables) sheetreplace


******************************************************
noi di _n(5) _dup(80) "=" _n "PHE discrenpencies significant" _n _dup(80) "="
********************************************************

cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Datasets"
use pipeline_clean_all_values_wide, clear

noi di "comparision between methods"
noi di  "Recall order is  Genefiner Mykrobe Typewriter"

gen agree  = (ValueA=="AAA"| ValueA=="PPP" )

*which set  have the most discrenpacies
noi di "site discrenpacies by set"
noi tab set agree, chi2 expected



* what is the minority method by site

gen minority= "GeneFinder" if agree!=1 &  valuet==valuez
replace minority= "Mykrobe" if agree!=1 &  valueg==valuet
replace minority= "Typewriter" if agree!=1 &  valueg==valuez

noi di "minority opinion by set, expected values included"
tab minority set, chi expected


*which set and type  have the most discrenpacies
noi di "site discrenpacies by set and type"
gen new = set + type
noi tab new minority, expected chi


noi di "method disagreeing with lab"

use anti_panel_all, clear
append using viru_panel_all, gen(new)




* reshape wide
drop _merge new valueall method
reshape long value, i (sample site set gold) j(method) string
replace value=upper(value)

drop if gold=="NA"
drop if gold=="I"
gen match=1 if gold==value
replace match=0 if match==.
tab set method if match!=1, chi

noi di "total number gold method results"

tab set method


******************************************************
noi di _n(5) _dup(80) "=" _n "Oxford validation/ training sets" _n _dup(80) "="
********************************************************

* t-tests of matches in two Oxford sets
ttest match if method=="zam" & set!="Collindale", by(set)
ttest match if method=="genefinder" & set!="Collindale", by(set)
ttest match if method=="typewriter" & set!="Collindale", by(set)


* sens and spec for each set individually


noi di "Collindale"
tempfile temp
	use anti_panel_all, clear
	assert site!=""
	assert _N==16548
	keep if set=="Collindale"
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
	gen set="Collindale"
save temp, replace
	
foreach k in valuetype valuez{
	use anti_panel_all, clear
	assert site!=""
	assert _N==16548
		keep if set=="Collindale"
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
	gen set="Collindale"
	append using temp
	save temp, replace
}

noi di "Oxford 491"
	foreach k in valuetype valuez valuegene{
	use anti_panel_all, clear
	assert site!=""
	assert _N==16548
		keep if set=="Oxford491"
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
	gen set="Oxford491"
	append using temp
	save temp, replace
}

noi di "Oxford 501"
	foreach k in valuetype valuez valuegene{
	use anti_panel_all, clear
	assert site!=""
	assert _N==16548
		keep if set=="Oxford501"
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
	gen set="Oxford501"
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
	cii  truePos[`i'] TP[`i']
	replace lsens = r(lb) if _n==`i'
	replace usens = r(ub) if _n==`i'
	cii truePos[`i'] FN[`i']
	replace lvme = r(lb) if _n==`i'
	replace uvme = r(ub) if _n==`i'
	}
	else{
	noi di site[`i'] " has no positive samples"
	}
	if trueNeg[`i']>0{
	cii trueNeg[`i'] TN[`i']
	replace lspec = r(lb) if _n==`i'
	replace uspec = r(ub) if _n==`i'
	cii trueNeg[`i'] FP[`i']
	replace lme = r(lb) if _n==`i'
	replace ume = r(ub) if _n==`i'
	}
	else{
	noi di site[`i'] " has no negative samples"
	}
	}

* save values for combining for table
noi list set method sens lsens usens spec lspec uspec
