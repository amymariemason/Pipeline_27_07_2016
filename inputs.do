* Import data from excel sheets
* based on code by Tim Peto

set li 130

cap log close
log using inputs.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Datasets"

************NEW DATA***************
*genefinder sets
noi di "Genefinder sets"
import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Genefinder.xlsx", sheet("Oxford-validation") firstrow allstring clear
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


import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Genefinder.xlsx", sheet("Oxford-derivation") firstrow allstring clear
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

import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Genefinder.xlsx", sheet("PHE-validation") firstrow allstring clear
gen method="genefinder"
gen set="Collindale"
drop CX
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
 import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\mykrobe_17_oct_2016.xlsx", sheet("mykrobe_17_oct_2016.tsv") firstrow clear

gen method="zam"
*gen set=""
*reshape to one line per site/sample
renpfix "" value
rename valuecomid sample
rename valuemethod method
rename valueremarks comments
reshape long value, i(sample method comments) j(site) string
drop comments
save pipeline_data_z, replace



* problem with extra samples; some in multiple sets with null results: drop any samples with duplicate results
noi di "duplicate results for same sample - one set indentical or null"

use pipeline_data_z, clear
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
import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Typewriter.xlsx", sheet("Colindale400_Michel_panel") cellrange(A2:CU400) firstrow allstring clear
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

import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Typewriter.xlsx", sheet("OxDerivation_Michel_panel") cellrange(A2:CU503) firstrow allstring clear

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

import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Typewriter.xlsx", sheet("OxValidation_Michel_panel") cellrange(A2:CU493) firstrow allstring clear
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
noi di "Lab gold standard"
import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Typewriter.xlsx", sheet("OxValidation_phenotype") cellrange(A2:O494) firstrow allstring clear
* picks up extra blank line - drop 
drop if A==""
gen set="Oxford491"
gen method="gold"
rename A sample
rename B comments
save pipeline_data_gs, replace

import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Typewriter.xlsx", sheet("OxDerivation_phenotype") cellrange(A2:M503) firstrow allstring clear

gen set="Oxford501"
gen method="gold"
rename A sample
append using pipeline_data_gs.dta, force
save pipeline_data_gs, replace


noi di "include extra toxin panal"
import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Genefinder.xlsx", sheet("PHE-phenotypes-and toxin PCRs") firstrow allstring clear
drop N
gen method="gold"
gen set="Collindale"
rename sample_id sample
append using pipeline_data_gs.dta, force

noi di _N " samples"

save pipeline_data_gs, replace

******************************************************************
* input coverage files 



* typewriter sets
noi di "typewriter coverage"

import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Typewriter_thresholds_10Feb15.xlsx", sheet("target list 1") firstrow clear
rename s twentythree_s
gen method="typewriter coverage"
rename BLASTDB sample
drop if sample ==""
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
reshape long value, i(sample method) j(site) string
save coverage_tw1, replace

import excel "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Inputs\Typewriter_thresholds_10Feb15.xlsx", sheet("target list 2") firstrow clear

gen method="typewriter coverage"
rename BLASTDB sample
drop if sample ==""
*reshape to one line per site/sample
renpfix "" value
rename valuesample sample
rename valuemethod method
reshape long value, i(sample method) j(site) string

append using coverage_tw1.dta, force
*clean sample names
gen break = strpos(sample, "_PHE")
gen shortsample = substr(sample, 1, break-1) if break!=0
replace shortsample=sample if break==0
drop sample
rename shortsample sample
save coverage_tw, replace

