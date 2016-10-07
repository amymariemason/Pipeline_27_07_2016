* Compare the three methods in phenotype predictions

set li 130

cap log close
log using analysis2.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

*********** create phenotype predictions by site
use anti_panel_all, clear

noi di" make bar chart of which methods differ on which antibiotic"
preserve
contract  site  value*
rename valueall all
drop value*

reshape wide _freq, i(site) j(all) string
rename _freq* *


* bar graph of predictions (antibiotics)
graph bar (asis)  rrr rrs rsr rss srs sss,  over(site, label(angle(90))) stack title("Combinations of Results") subtitle("Results given as Genefinder Mykrobe Typewriter ") legend( label(1 "rrr") label( 2 "rrs" ) label (3 "rsr") label (4 "rss" ) label (5 "srs") label( 6 "sss" )) bar(1, color(gs8)) bar(2, color(cyan)) bar(3, color(red)) bar(4, color(blue)) bar(5, color(yellow)) bar(6, color(gs12))
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\phenotype_disagreements.tif", as(tif) width(2550) replace

restore
*************************************************************

noi di "compare to gold standard"

* in table
noi tab gold valueall

* in graph
preserve
contract  site  gold valueall

reshape wide _freq, i(site gold) j(valueall) string
rename _freq* *


gen graphx = site + " " + gold
drop if !inlist(gold, "R", "S")
sort gold site
 gen num =_n
 replace num=num+1 if gold=="S"
labmask num, values(site)
gen goldlabel = "Gold standard " + gold

* bar graph of predictions (antibiotics)
graph bar (asis) rrr rrs rsr rss srs sss if inlist(golds,"R", "S"),graphregion(color(white))    over(site, label(angle(90) labsize(tiny)) )  over(goldlabel, label(labsize(small)) ) stack  title("Prediction combinations") subtitle("Results given as Genefinder Mykrobe Typewriter") legend(rows(2) label(1 "rrr") label( 2 "rrs" ) label (3 "rsr") label (4 "rss" ) label (5 "srs") label( 6 "sss" )) bar(1, color(gs8)) bar(2, color(cyan)) bar(3, color(red)) bar(4, color(blue)) bar(5, color(yellow)) bar(6, color(gs12)) ylabel(0(500)1000 1400)
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\phenotype_disagreements_withgold.tif", as(tif) width(2550) replace

restore

* sensitivity/ specificity
noi di "find sensitivity and specifity of each method/antibiotic"

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

preserve
keep site  sens spec lsens usens lspec uspec
rename * tw_*
rename tw_site site
save typewriter, replace
restore

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
graph save Graph "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_sens.gph", replace
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\tw_sens.tif", as(tif) replace

