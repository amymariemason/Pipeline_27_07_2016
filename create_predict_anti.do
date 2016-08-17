
* create resistance prediction panal for all methods; add goldstandard

set li 130

cap log close
log using panal_anti.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

************************
* add antibiotic relevance to sites

******************************************************
noi di _n(5) _dup(80) "=" _n " 1 create antibiotic resistance prediction panels" _n _dup(80) "="
**************************
************************
use pipeline_clean_all_values_long, clear

gen antibiotic =""
replace antibiotic = "ciprofloxacin" if inlist(site, "gyra","gryb", "grla", "grlb")
replace antibiotic = "erythromycinclindamycin" if inlist(site, "erma", "ermb", "ermc", "ermt", "ermy")
replace antibiotic = "erythromycin" if inlist(site, "msra", "mphc")
replace antibiotic = "fusidicacid" if inlist(site, "fusb", "fusc", "fusa", "far")
replace antibiotic = "gentamicin" if inlist(site,"aac6aph2", "aph2ic")
replace antibiotic = "methicillinpenicillin" if inlist(site,"meca", "mecc")
replace antibiotic = "mupirocin" if inlist(site, "mupa", "mupb")
replace antibiotic = "penicillin" if site=="blaz"
replace antibiotic = "rifampicin" if site=="rpob"
replace antibiotic = "tetracycline" if inlist(site, "tetk", "tetl", "tetm", "teto")
replace antibiotic = "trimethoprim" if inlist(site, "dfra", "dfrc", "dfrb", "dfrb","dfrg", "dfrk")
replace antibiotic = "vancomycin" if inlist(site, "vana", "vanb", "vanc")

******* keep just antibiotics

keep if antibiotic!=""
compress
gen marker = (value=="P")
drop site set type

* collapse by antibiotic, so that one record per antibiotic/sample/method
bysort sample method antibiotic: egen maxmarker= max(marker)
duplicates drop sample method antibiotic maxmarker, force
drop value marker


* now create wide record, one column per antibiotic
reshape wide maxmarker, i( sample method) j(antibiotic) string


* clean up multiple indicators
rename maxmarker* *
replace erythromycin = 1 if (erythromycinclindamycin==1)
rename erythromycinclindamycin clindamycin
gen methicillin =  methicillinpenicillin
replace penicillin = 1 if  methicillinpenicillin==1
drop methicillinpenicillin


foreach k in  ciprofloxacin erythromycin clindamycin fusidicacid gentamicin mupirocin penicillin rifampicin tetracycline trimethoprim vancomycin methicillin{
		noi di "`k'"
		gen predict`k'= "r" if `k'==1
		replace predict`k' = "s" if `k'==0
		drop `k'
		rename predict`k' `k'
		noi tab `k', m
		}

noi display "save panal"
save anti_prediction, replace

preserve

rename * value*
rename valuesample sample
rename valuemethod method
reshape long value, i(sample method) j(site) string
noi display "basic antibiotic prediction counts"
noi bysort method: tab site value 
noi bysort site: tab method value
noi display "counts differing between method: gene, type, mykrobe"

reshape wide value, i(sample site) j(method) string
gen valueall = valueg + valuet + valuez
noi tab valuea, m sort

noi tab site valuea, m 
save anti_prediction_long, replace

restore
**********************************

******************************************************
noi di _n(5) _dup(80) "=" _n " 1 add gold standard" _n _dup(80) "="
**************************
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets
use pipeline_gold_clean_long, clear
drop if inlist(upper(site), "SEA", "SEB", "SEC", "SED", "SEE", "SEG", "SEH")
drop if inlist(upper(site), "SEI", "SEJ", "SELR")
drop if inlist(upper(site),  "SEP", "SEU", "ETA", "ETB", "ETD", "TSST1", "PVL")
drop if inlist(upper(site), "MECA", "MECC")
noi merge 1:1 sample site using anti_prediction_long, update
assert _merge==3

noi di "predictive is gene, type, mykrobe"
noi tab gold valuea

noi bysort gold: tab site valuea

save anti_panel_all, replace
cd E:\users\amy.mason\Pipeline_27_07_2016\

