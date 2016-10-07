set li 130

cap log close
log using mykrobe_prob.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

* create list of mykrobe disagreements for Phelim

use anti_panel_all, clear

keep if valueall=="rsr" | valueall == "srs"
noi di _N "sample- antibiotic predictions where mykrobe disgrees with typewriter and genefinder"

noi bysort valuezam: tab site gold 
noi tab set
 export excel using "E:\users\amy.mason\Pipeline_27_07_2016\Excel_Outputs\Mykrobe_disagreements.xls", replace firstrow(variables)

 **********
 * why is ciprofloxin not more of a problem - because most of the samples that are "PAP" for gyla are "PPP" for grla (207/212)
 
use pipeline_clean_all_values_wide, clear
 keep if ValueAll=="APA" | ValueAll=="PAP"

  export excel using "E:\users\amy.mason\Pipeline_27_07_2016\Excel_Outputs\Mykrobe_site_disagreements.xls", replace firstrow(variables)

 
 noi di  _N "site-results where mykrobe disgrees with typewriter and genefinder"
noi bysort valuez: tab site

noi bysort site: tab sample ValueAll

keep if site=="gyra"
keep sample
save my_prob, replace

use pipeline_clean_all_values_wide, clear
merge m:1 sample using "my_prob", update
keep if _merge==3

keep if inlist(site, "gyra","gryb", "grla", "grlb")

noi tab site ValueAll
noi di "207/ 212 sites where myrkobe didn't identify gyra, did find grla"
*************************************************************