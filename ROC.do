************************************************
* ROC.DO 
************************************************
************************************************
* Creates ROC curves to analysis optimum thresholds from the Typewriter data

*Inputs: coverage_tw (from INPUTS.DO),  pipeline_gold_clean_wide (from clean_pheno.do)
* Outputs : ROC curves (see folder ROC inside GRAPHS)


* Written by: Amy Mason


set li 130

cap log close
log using ROC.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets



* Create a binary gold standard set; (swap R/S to 0/1)

use pipeline_gold_clean_wide, clear
rename * gold*
rename (goldsample goldmethod goldset) (sample method set)
reshape long gold, i(sample) j(anti) string
gen value= 1 if gold=="R"
replace value = 0 if gold=="S"
drop gold
reshape wide value, i(sample) j(anti) string
rename value* gold*
save gold, replace

* merge lab data and coverage data
use coverage_tw, clear
reshape wide
drop method
merge 1:1 sample using gold, update
* not in using if dropped samples
drop if _merge==1 
assert _merge==3
drop _merge
rename value* *

save temp, replace

*******************************************************************************
* Create ROC Curves and identify optimal cut-off for each antibiotic
*************************************************************************
*************
* clindamycin
**************
* create ROC graph
use temp, clear

egen clindamax = rowmax(ermA_1 ermA_2 ermB ermC ermY ermT)
roctab goldClinda clindamax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\clinda.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldClinda clindamax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of ermA_1 ermA_2 ermB ermC ermY ermT in predicting clindamycin resistance"
noi list clindamax sens spec 

*********************************************************
* Ciprofloxacin
*********************************************************
* complicated as combination prediction of chromosonal mutations; unsure how to predict

*********************************************************
* Erythromycin
********************************************************* 
use temp, clear
egen erythmax= rowmax(ermA_1 ermA_2 ermB ermC ermY ermT msrA mphC)
roctab goldEryth erythmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\Erythromycin.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldEryth erythmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of ermA_1 ermA_2 ermB ermC ermY ermT msrA mphC in predicting Erythromycin resistance"
noi list erythmax sens spec 

 
 ****************************************************************
 *gentamicin
 ****************************************************************
use temp, clear
egen gentamax= rowmax(aac6aph2 aph2Ic)
roctab goldGenta gentamax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\gentamicin.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldGenta gentamax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of aac6aph2 aph2Ic in predicting gentamicin resistance"
noi list gentamax sens spec 
 
***************************************************************
*methicillin
****************************************************************
use temp, clear
egen methmax= rowmax(mecA mecC)
roctab goldMeth methmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\methicillin.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldMeth methmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of aac6aph2 aph2Ic in predicting methicillin resistance"
noi list methmax sens spec 
 
 *******************************************************************
*Penicillin
************************************************************************
* Problem: mecA and MecC have cutoffs of 80, blaz of 30 -  rescale?
* first: do as previous, just find max
use temp, clear
egen penmax= rowmax(mecA mecC blaZ)
roctab goldPen penmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\Penicillin.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldPen penmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of mecA mecC blaz in predicting Penicillin resistance"
noi list mecA mecC blaZ penmax sens spec 

*************************************** rescale attempt
*mecA and MecC have cutoffs of 80, blaz of 30 -  rescale?
use temp, clear
gen blaztemp = blaZ/30*80
replace blazt = 100 if blazt>100
egen penmax= rowmax(mecA mecC blazt)
roctab goldPen penmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\Penicillin_rescale.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldPen penmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of mecA mecC  and blaz/30*80 (with ceiling at 100) in predicting Penicillin resistance"
noi list mecA mecC blaZ penmax sens spec 

************* blaZ only 
use temp, clear
egen penmax= rowmax(blaZ)
roctab goldPen penmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\Penicillin_blaZ_only.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldPen penmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of blaZ in predicting Penicillin resistance"
noi list mecA mecC blaZ penmax sens spec 

 ********************************************
* Mupirocin
*********************************************
use temp, clear
egen mupmax= rowmax( mupA  mupB)
roctab goldMup mupmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\Mupirocin.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldMup mupmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert mupmax == mupmax[_n-1] if _n>1
keep if _n==1

noi di "optimum cut-off for maximum of mupA  mupB in predicting Mupirocin resistance"
noi list mupmax sens spec 

**************************************************************
*tetracycline
**************************************************

use temp, clear
egen tetramax= rowmax( tetK tetL tetM tetO)
roctab goldTetra tetramax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\tetracycline.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldTetra tetramax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of tetK tetL tetM tetO in predicting tetracycline resistance"
noi list tetramax sens spec 
 
************************************************
*vancomycin
*************************************************
* doesn't work as all lab results sensitive

***************************************************************
* *fusidicacid
**************************************************************
* this has a combination of chromosonal and aquired resistance. Need to get rid of the samples positive for mutations at fusa and analysis in their absense

use pipeline_clean_all_values_long, clear
keep if method=="typewriter"
keep if site=="fusa"
merge 1:1 sample using temp, update
drop if value=="P"
egen fusmax= rowmax( fus_B fus_C)
roctab goldFusi fusmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\fusidicacid.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldFusi fusmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert _n==1

noi di "optimum cut-off for maximum of fus_B fus_C in absense of positive fusA mutation in predicting fusidicacid resistance"
noi list fusmax sens spec 
 
 **************************************************************************
  *rifampicin
  ****************************************************************
  *only chromosonal resistance,
 
 ****************************************************************
 *trimethoprim
 ********************************************************************
use pipeline_clean_all_values_long, clear
keep if method=="typewriter"
keep if site=="dfrb"
merge 1:1 sample using temp, update
drop if value=="P"
egen trimax= rowmax( dfr_A dfr_C dfr_D dfr_G dfr_K)
roctab goldTri trimax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\trimethoprim.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldTri trimax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
assert trimax == trimax[_n-1] if _n>1
keep if _n==1

noi di "optimum cut-off for maximum of dfr_A dfr_C dfr_D dfr_G dfr_K in absense of positive dfr_B mutation in predicting fusidicacid resistance"
noi list trimax sens spec 
 
 
 *******************************************
  ****************************************************************

 * Virulence panel

  ****************************************************************
*******************************************


********************************************************** 
 *pvl
 *********************************************************
 
 
 use temp, clear
egen lukmax= rowmax( lukPVF lukPVS)
roctab goldPVL lukmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\pvl.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldPVL lukmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
keep if _n==1

noi di "optimum cut-off for maximum of lukPVF lukPVS in predicting PVL virulence"
noi list lukmax sens spec 
 
 
  ****************************************************************
* SEA 
 ****************************************************************
 
  use temp, clear
egen seamax= rowmax(sea_1 sea_2)
roctab goldSEA seamax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\SEA.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldSEA seamax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
keep if _n==1

noi di "optimum cut-off for maximum of sea_1 sea_2 in predicting SEA virulence"
noi list seamax sens spec 
 
 
 
 
 ****************************************************************
 * SEB -SEP
  ****************************************************************
  *single sites, so do as a loop
  
foreach x in seb sec see seg sei sej eta etb etd tsst1 meca mecc{ 
  use temp, clear
  rename *,lower
roctab gold`x' `x', graph
cd E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC
graph export `x'.tif, as(tif) replace
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec gold`x' `x', sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
keep if _n==1

noi di "optimum cut-off for maximum of " `x' " in predicting " `x' " virulence"
noi list `x' sens spec 
 
 }
 
  ****************************************************************
 * SEH 
  ****************************************************************
  use temp, clear
egen sehmax= rowmax(seh_1 seh_2)
roctab goldSEH sehmax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\SEH.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldSEH sehmax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
keep if _n==1

noi di "optimum cut-off for maximum of seh_1 seh_2 in predicting SEH virulence"
noi list sehmax sens spec 
 
 
 ****************************************************************
 * SEU 
 ****************************************************************
  use temp, clear
egen seumax= rowmax(seu_1 seu_2)
roctab goldSEU seumax, graph
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\ROC\SEU.tif", as(tif) replace


*use senspec to identify optimal cut-off via Youden/ distance arguments
senspec goldSEU seumax, sensitivity(sens) specificity(spec)

gen youden= sens-(1-spec)
egen youdenmax= max(youden)
gen count = (youden==youdenmax)
keep if count ==1
keep if _n==1

noi di "optimum cut-off for maximum of seu_1 seu_2 in predicting SEU virulence"
noi list seumax sens spec 
 
 
