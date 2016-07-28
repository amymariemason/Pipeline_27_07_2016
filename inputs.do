* Import data from excel sheets
* based on code by Tim Peto

set li 130

cap log close
log using inputs.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

************NEW DATA***************
*genefinder sets
noi di "Genefinder sets"
import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Genefinder_pipeline_conparison_Dec-2015.xlsx", sheet("Oxford-validation") firstrow allstring clear
gen method="genefinder"
gen set="Oxford491"
rename Sample_id sample
rename remarks comments
* there is a chromo site called 23s/ s, only genefinder reports on it: drop it
drop s
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
save pipeline_data_gf, replace


import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Genefinder_pipeline_conparison_Dec-2015.xlsx", sheet("Oxford-derivation") firstrow allstring clear
gen method="genefinder"
gen set="Oxford501"
rename Sample_id sample
rename remarks comments
* there is a chromo site called 23s/ s, only genefinder reports on it: drop it
drop s
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
append using pipeline_data_gf.dta, force
save pipeline_data_gf, replace

import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Genefinder_pipeline_conparison_Dec-2015.xlsx", sheet("PHE-validation") firstrow allstring clear
gen method="genefinder"
gen set="Collindale"
rename Sample_id sample
rename remarks comments
* there is a chromo site called 23s/ s, only genefinder reports on it: drop it
drop s
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
append using pipeline_data_gf.dta, force
noi di _N " results"
bysort site: gen tab=1 if _n==1
summ tab 
noi di "on " r(sum) " sites"
drop tab
bysort sample: gen tab=1 if _n==1
summ tab
noi di "from " r(sum) " samples"
drop tab
save pipeline_data_gf, replace



* Mykroke sets
noi di "Mykrobe sets"
import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\mykrobe_d425f57e5c526aec54eed1a51a875b82483998c3.xlsx", sheet("Oxford 501 derivation") cellrange(A2:CJ544) firstrow allstring  clear
gen method="zam"
gen set="Oxford501"
rename A sample
rename B comments
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
save pipeline_data_z, replace

import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\mykrobe_d425f57e5c526aec54eed1a51a875b82483998c3.xlsx", sheet("Oxford 491 validation") cellrange(A2:CJ544) firstrow allstring  clear
gen method="zam"
gen set="Oxford491"
rename A sample
rename B comments
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
append using pipeline_data_z.dta, force
save pipeline_data_z, replace

import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\mykrobe_d425f57e5c526aec54eed1a51a875b82483998c3.xlsx", sheet("Colindale 400") cellrange(A2:CJ544) firstrow allstring clear
gen method="zam"
gen set="Collindale"
rename A sample
rename B comments
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
append using pipeline_data_z.dta, force
save pipeline_data_z, replace


* problem with extra samples; some in multiple sets with null results: drop any samples with duplicate results
noi di "duplicate results for same sample - one set indentical or null"

use pipeline_data_z, clear
gsort sample site value -set
by sample site value: gen extra = 1 if _n>1 
replace extra=1 if inlist(sample, "C00011656", "C00011655") & comments==""
replace extra=1  if sample=="C00007108" & strpos(set, "dale")
summ extra if extra==1
noi di r(N) " duplicate results dropped"
noi drop if extra==1
drop extra 
noi di _N " results"
bysort site: gen tab=1 if _n==1
summ tab 
noi di "on " r(sum) " sites"
drop tab
bysort sample: gen tab=1 if _n==1
summ tab
noi di "from " r(sum) " samples"
drop tab
save pipeline_data_z2, replace

* typewriter sets
noi di "typewriter"
import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Pipeline_comparison_Typewriter_results_Michel_panel_23Jan15(1).xlsx", sheet("Colindale400_Michel_panel") cellrange(A2:CU400) firstrow allstring clear
* picks up extra blank line from excel  - drop
drop if Sample ==""
gen set="Collindale"
gen method="typewriter"
rename SampleID sample
rename B comments
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
save pipeline_data_tw, replace

import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Pipeline_comparison_Typewriter_results_Michel_panel_23Jan15(1).xlsx", sheet("OxDerivation_Michel_panel") cellrange(A2:CU503) firstrow allstring clear

gen set="Oxford501"
gen method="typewriter"
rename SampleID sample
rename B comments
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
append using pipeline_data_tw.dta, force
save pipeline_data_tw, replace

import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Pipeline_comparison_Typewriter_results_Michel_panel_23Jan15(1).xlsx", sheet("OxValidation_Michel_panel") cellrange(A2:CU493) firstrow allstring clear
gen set="Oxford491"
gen method="typewriter"
rename SampleID sample
rename B comments
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
append using pipeline_data_tw.dta, force
noi di _N " results"
bysort site: gen tab=1 if _n==1
summ tab 
noi di "on " r(sum) " sites"
drop tab
bysort sample: gen tab=1 if _n==1
summ tab
noi di "from " r(sum) " samples"
drop tab
save pipeline_data_tw, replace

***************** import gold standard 

******* WHY IS THE METHOD HERE NOT GOLD STANDARD?


noi di "Lab gold standard"
import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Pipeline_comparison_Typewriter_results_Michel_panel_23Jan15(1).xlsx", sheet("OxValidation_phenotype") cellrange(A2:O494) firstrow allstring clear
* picks up extra blank line - drop 
drop if A==""
gen set="Oxford491"
gen method="gold"
rename A sample
rename B comments
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
rename valuecomments comments
reshape long value, i(sample set method comments) j(site) string
save pipeline_data_gs, replace

import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Pipeline_comparison_Typewriter_results_Michel_panel_23Jan15(1).xlsx", sheet("OxDerivation_phenotype") cellrange(A2:M503) firstrow allstring clear

gen set="Oxford501"
gen method="gold"
rename A sample
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
reshape long value, i(sample set method) j(site) string
append using pipeline_data_gs.dta, force
save pipeline_data_gs, replace

import excel "E:\users\amy.mason\Pipeline_27_07_2016\Inputs\Genefinder_pipeline_conparison_Dec-2015.xlsx", sheet("PHE-phenotypes-PVL") firstrow allstring clear
drop O
drop PVLbyPCR
gen method="gold"
gen set="Collindale"
rename sample_id sample
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
rename valueset set
reshape long value, i(sample set method) j(site) string
append using pipeline_data_gs.dta, force

noi di _N " results"
bysort site: gen tab=1 if _n==1
summ tab 
noi di "on " r(sum) " sites"
drop tab
bysort sample: gen tab=1 if _n==1
summ tab
noi di "from " r(sum) " samples"
drop tab

save pipeline_data_gs, replace



****** combine the 4 sets
noi di "combine"
use pipeline_data_gs, clear
append using pipeline_data_gf.dta, force
append using pipeline_data_tw.dta, force
append using pipeline_data_z2.dta, force
replace site = lower(site)
noi di _N " results"
bysort site: gen tab=1 if _n==1
summ tab 
noi di "on " r(sum) " sites"
drop tab
bysort sample: gen tab=1 if _n==1
summ tab
noi di "from " r(sum) " samples"
drop tab


cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets
save pipeline_all, replace
cd E:\users\amy.mason\Pipeline_27_07_2016\

