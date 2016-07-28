***************
* consider the sensitivity and specificity of all three methods combined against lab results
*run analysis_phenotypes first
****************

* upload lab data

cd E:\users\amy.mason\Pipeline\Amy_work

use "E:\users\amy.mason\Pipeline\Amy_work\lab_compare_all_wide.dta", clear

cd E:\users\amy.mason\Pipeline_27_07_2016\


* get wetlab counts
reshape long

* add type to discriminate graphs at later stage
gen type =""
replace type ="Phenotype" if inlist(site, "ciprofloxacin", "fusidicacid", "clindamycin", "erythromycin")
replace type ="Phenotype" if inlist(site, "gentamicin", "mupirocin", "penicillin","rifampicin")
replace type ="Phenotype" if inlist(site,"tetracycline", "trimethoprim", "vancomycin", "methicillin")

keep if type!=""
*keep if type==""

gen anti=(strlen(site)>5)

contract  anti site gold gene typew zam
gen result = gene + typew + zam 

* note order is GENEFINDER THEN TYPEWRITER THEN MYKROBE
drop gene typew zam
replace result = "missing" if result==""
replace gold = "missing" if gold==""


reshape wide _freq, i(site gold) j(result) string
rename _freq* *

egen combor = rowtotal(RRR RSR)
replace combor=0 if combor==.
egen combos= rowtotal(SSR SSS)
replace combos=0 if combos==.

 drop R* S*
 
 keep if inlist(gold, "R", "S")
 
 
 reshape wide combo*, i(site) j(gold) string
 
 rename combo* *

 *** generate sensitivity and specificity
 ***** NOTE: [ called GOLD] is format of these variables

 replace sS = 0 if sS==.
 replace rS = 0 if rS==.
 replace sR = 0 if sR==.
 replace rR = 0 if rR==.
   

for any lsens usens lspec uspec lme ume lvme uvme:gen X=.

gen TP= rR
gen FP = rS
gen TN = sS
gen FN = sR

gen POS = TP+FN
gen NEG = TN+ FP
gen specificity =TN/(TN+FP)
gen sensitivity = TP/(TP+FN)
gen MajorErrorRate = FP / (FP + TN)
gen VeryMajorErrorRate = FN/(TP+FN)
gen Accuracy = (TN + TP)/(TN+TP+FN+FP)

local max=_N-1
forvalues i=1(1)`max'{
di site[`i']
cii prop NEG[`i'] TN[`i']
replace lspec = r(lb) if _n==`i'
replace uspec = r(ub) if _n==`i'
cii prop POS[`i'] TP[`i']
replace lsens = r(lb) if _n==`i'
replace usens = r(ub) if _n==`i'
cii prop NEG[`i'] FP[`i']
replace lme = r(lb) if _n==`i'
replace ume = r(ub) if _n==`i'
cii prop POS[`i'] FN[`i']
replace lvme = r(lb) if _n==`i'
replace uvme = r(ub) if _n==`i'
}

* vancomycin has no
di site[_N]
cii prop NEG[_N] TN[_N]
replace lspec = r(lb) if _n==_N
replace uspec = r(ub) if _n==_N

cii prop NEG[_N] FP[_N]
replace lme = r(lb) if _n==_N
replace ume = r(ub) if _n==_N


sort site
 gen num =_n 

 order site  Major lme ume Very lvme uvme sS sR rS rR
 format  MajorErrorRate lme ume VeryMajorErrorRate lvme uvme lsens usens lspec uspec sensitivity specificity Accuracy num %4.3f



labmask num, values(site)
preserve
rename * maj_*
rename maj_site site
save majority, replace
restore

local max= _N
#delimit ;
twoway rcap usens lsens num || scatter sens num,  xlabel(1(1)`max', valuelabel angle(90))
legend(off) title("Majority Rule vs. Lab") 
subtitle("Sensitivity in Phenotype prediction");
#delimit cr
graph save Graph "E:\users\amy.mason\Pipeline_27_07_2016\majority_sens_all.gph", replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\majority_sens_all.tif", as(tif) replace


local max= _N
#delimit ;
twoway rcap lspec uspec num || scatter spec num,   xlabel(1(1)`max', valuelabel angle(90)) 
legend(off) title("Majority Rule vs. Lab") 
subtitle("Specifitity in Phenotype prediction")
xtitle("Gene site");
#delimit cr

graph save Graph "E:\users\amy.mason\Pipeline_27_07_2016\majority_spec_all.gph", replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\majority_spec_all.tif", as(tif) replace




collapse (sum) sS sR rS rR, by(anti)

for any lsens usens lspec uspec lme ume lvme uvme:gen X=.


gen cPOS = sS+sR
gen cNEG = rR +rS
gen specificity =sS/cPOS
gen sensitivity = rR/cNEG
gen truNeg= rR+sR
gen truPos= sS+ rS
gen MajorErrorRate = rS/truPos
gen VeryMajorErrorRate = sR/truNeg
gen Agreement = (sS+rR)/(truNeg+truPos)

local max=_N
forvalues i=1(1)`max'{
cii  cPOS[`i'] sS[`i']
replace lspec = r(lb) if _n==`i'
replace uspec = r(ub) if _n==`i'
cii  cNEG[`i'] rR[`i']
replace lsens = r(lb) if _n==`i'
replace usens = r(ub) if _n==`i'
cii  truPos[`i'] rS[`i']
replace lme = r(lb) if _n==`i'
replace ume = r(ub) if _n==`i'
cii  truNeg[`i'] sR[`i']
replace lvme = r(lb) if _n==`i'
replace uvme = r(ub) if _n==`i'
}


collapse (sum) sS sR rS rR

for any lsens usens lspec uspec lme ume lvme uvme:gen X=.


gen cPOS = sS+sR
gen cNEG = rR +rS
gen specificity =sS/cPOS
gen sensitivity = rR/cNEG
gen truNeg= rR+sR
gen truPos= sS+ rS
gen MajorErrorRate = rS/truPos
gen VeryMajorErrorRate = sR/truNeg
gen Agreement = (sS+rR)/(truNeg+truPos)

local max=_N
forvalues i=1(1)`max'{
cii  cPOS[`i'] sS[`i']
replace lspec = r(lb) if _n==`i'
replace uspec = r(ub) if _n==`i'
cii  cNEG[`i'] rR[`i']
replace lsens = r(lb) if _n==`i'
replace usens = r(ub) if _n==`i'
cii  truPos[`i'] rS[`i']
replace lme = r(lb) if _n==`i'
replace ume = r(ub) if _n==`i'
cii  truNeg[`i'] sR[`i']
replace lvme = r(lb) if _n==`i'
replace uvme = r(ub) if _n==`i'
}
******** create combined graphs ****

use majority, clear
merge 1:1 site using typewriter
drop if _merge!=3
drop _merge
merge 1:1 site using genefinder
drop if _merge!=3
drop _merge
merge 1:1 site using mykrobe
drop if _merge!=3
drop _merge

* get into right shape for graph

reshape long maj_ type_ gene_ myk_, i(site) j(new) string
rename maj_ meth_maj
rename type_ meth_type
rename gene_ meth_gene
rename myk_ meth_myk
reshape long meth_, i(site new) j(method) string
rename meth_ value
reshape wide value, i(site method) j(new) string
rename value* *

* sneaky creating gaps on graphs
expand 4 if method =="maj", gen(dups)
replace method ="blank" if dups==1
replace sens=. if method=="blank"
replace lsens=. if method=="blank"
replace usens=. if method=="blank"
replace uspec=. if method=="blank"
replace lspec=. if method=="blank"
replace spec=. if method=="blank"


encode(method), gen(order)
replace order = 8 if order==2
sort site order
gen newcount = _n
local max =_N
labmask newcount, values(site)

#delimit ;
twoway (rcapsym usens lsens newcount, s(i) ylabel(,format(%3.2f))) (scatter sens newcount if method=="blank") (scatter  sens newcount if method=="gene", mcolor(orange))  (scatter  sensitivity newcount if method =="type", mcolor(red)) 
(scatter  sens newcount if method=="myk", mcolor(blue))  (scatter sens newcount if method=="maj", msymbol(D) mcolor(black)),  xlabel(5(7)84, valuelabel angle(90))
legend(order(3 4 5 6) lab(3 "Genefinder") lab(4 "Typewriter") lab(5 "Mykrobe") lab(6 "Majority Rule"))
subtitle("Sensitivity in Phenotype prediction") graphregion(fcolor(white));
#delimit cr

graph save Graph "E:\users\amy.mason\Pipeline_27_07_2016\Sens_all.gph", replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Sens_all.tif", as(tif) replace


#delimit ;
twoway (rcapsym uspec lspec newcount,s(i)  ylabel(,format(%3.2f))) (scatter sens newcount if method=="blank") (scatter  spec newcount if method=="gene", mcolor(orange))  (scatter  spec newcount if method =="type", mcolor(red)) 
(scatter  spec newcount if method=="myk", mcolor(blue))  (scatter spec newcount if method=="maj", msymbol(D) mcolor(black)),  xlabel(5(7)84, valuelabel angle(90))
legend(order(3 4 5 6) lab(3 "Genefinder") lab(4 "Typewriter") lab(5 "Mykrobe") lab(6 "Majority Rule"))
subtitle("Specificity in Phenotype prediction") graphregion(fcolor(white));
#delimit cr

graph save Graph "E:\users\amy.mason\Pipeline_27_07_2016\Spec_all.gph", replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Spec_all.tif", as(tif) replace


