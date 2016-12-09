************************************************
* CREATE_PREDICT_VIRU.DO 
************************************************

* Takes the  site predictions from the three methods, and creates virulence predictions

*Inputs:  pipeline_clean_all_values_long (one record per sample, per method, per site) (from clean_predict.do) and pipeline_gold_clean_long (from clean_pheno.do)
* Outputs : virulence_prediction_long  (one record per sample per virulence factor)  viru_panel_all (one record per sample per virulence factor, with the three predictions and the gold standard result if known)


* Written by: Amy Mason





set li 130

cap log close
log using panal_viru.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

************************
* VIRULENCE PREDICTION PANEL
******************************************************
noi di _n(5) _dup(80) "=" _n " 1 create virulence prediction panels" _n _dup(80) "="
**************************
* for each virulence factor in our panel, and for each sample we want to know which methods predict it present or not and what the final lab result was.
************************
use pipeline_clean_all_values_long, clear

* keep only the sites that relevant to virulence

no di "drop sites with no virulence test"
keep if strpos(site, "tsst") | strpos(site, "et") | strpos(site, "se") | strpos(site, "luk") | strpos(site, "mec")
drop if site=="lukm" | site=="lukmf"
drop if strpos(site, "tet") 
compress
noi tab site, sort

* create presence marker
gen marker = (value=="P")
drop set type value


* collapse by site, so that one record per antibiotic/sample/method
duplicates drop sample method site, force
bysort sample method site: assert _n==1



* now create wide record, one column per antibiotic
reshape wide marker, i(sample method) j(site) string
rename marker* *

* clean up multiple indicators
gen pvl = 1 if  lukpvf==1 & lukpvs ==1
replace pvl =0 if lukpvf==0 | lukpvs ==0
drop luk*

*generate prediction variable; use a/p as resistance/ sensitive doesn't make sense in this context
foreach k in   eta etb etd meca mecc sea seb sec sed see seg seh sei sej selr sep seu tsst1 pvl{
		noi di "`k'"
		gen predict`k'= "p" if `k'==1
		replace predict`k' = "a" if `k'==0
		drop `k'
		rename predict`k' `k'
		noi tab `k', m
		}

		* save
noi display "save panal"

save virulence_prediction, replace

* reshape long and save
rename * value*
rename valuesample sample
rename valuemethod method
reshape long value, i(sample method) j(site) string
noi display "basic virulence prediction counts"
noi bysort method: tab site value 
noi bysort site: tab method value
noi display "counts differing between method: gene, type, mykrobe"

reshape wide value, i(sample site) j(method) string
gen valueall = valueg + valuez + valuet
noi tab valuea, m sort

noi tab site valuea, m 
save virulence_prediction_long, replace


******************************************************

*****************************************************
noi di _n(5) _dup(80) "=" _n " 1 add gold standard" _n _dup(80) "="
****************************************************
* add in the gold standard results from clean_pheno; note not all values have results

cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets
use pipeline_gold_clean_long, clear
drop if inlist(site, "ciprofloxacin","erythromycin","clindamycin","fusidicacid","gentamicin")
drop if inlist(site,"mupirocin","penicillin","rifampicin","tetracycline")
drop if inlist(site,"trimethoprim","vancomycin","methicillin")
		
* swap to A/P as resistance/ sensitive doesn't make sense in this context		
replace gold = "P" if gold=="R"
replace gold ="A" if gold=="S"
		
noi merge 1:1 sample site using virulence_prediction_long, update
assert _merge==3

noi di "predictive is gene,  mykrobe, type,"
noi tab gold valuea

noi bysort gold: tab site valuea
save viru_panel_all, replace
cd E:\users\amy.mason\Pipeline_27_07_2016\