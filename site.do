************************************************
*SITE.DO 
************************************************

* It was suggested that the discrepent results might be due to the results being close to threshold levels; This analysis picks out the coverage values of the discrepent results
* and looks at the coverage levels of those results. There is no evidence to suggest that many discrepent results are due to threshold coverage  as < 10% of them were within
* the 75-85 range (25-35 for blaz).

*Inputs:pipeline_clean_all_values_wide (from clean_predict.do), coverage_tw (from input.do) , gene_coverate_zam (from input.do)
* Outputs : none

* Written by: Amy Mason

************************************************************************

************************************************************************

set li 130

cap log close
log using site.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

use pipeline_clean_all_values_wide, clear

**************************
* find the sites where the values do not match between all three methods
noi di  "Recall order is  Mykrobe Mykrobe Genefinder"

gen agree  = (ValueA=="AAA"| ValueA=="PPP" )
drop if agree==1
noi di "disagreements: sites"
noi tab site, sort
noi di "disagreements: samples"
noi tab sample, sort
noi di "disagreements: set"
noi tab set, sort
noi tab site set

save discrepancies, replace
local samplesites = _N
di _N

**********************************
* Compare to typewriter values : cleaning and matching
********************************************************************
noi di _n(5) _dup(80) "=" _n " Typewriter values" _n _dup(80) "="
********************************************************************

use coverage_tw, clear

* as with the original a/p data, there is some cleaning up that needs to be done to match site names
* this is based on the cleaning code in clean_predict.do

* need to clean site name values 
**** match up things I think are the same but with slightly different names
replace site = lower(site)
replace site = "ant9ib" if site=="ant9_ibspc"
replace site = "ant9ia" if site=="ant9_iaspc"
replace site= "aac6aph2" if strpos(site, "aac")

replace site = subinstr(site, "acr", "arc",.) if strpos(site, "acr")
replace site = "mlst" if site=="st"
replace site = "ccrca"  if site == "ccrc_a"
replace site = "ccrcb"  if site == "ccrc_b"
replace site = "ccrcc"  if site == "ccrc_c"
replace site =  subinstr(site, "_", "",.)  if inlist(site, "dfr_a", "dfr_c", "dfr_d", "dfr_g", "dfr_k")
replace site =  subinstr(site, "_", "",.)  if inlist(site, "fus_b", "fus_c")
replace site =  subinstr(site, "_", "",.)  if strpos(site, "qac")
replace site =  subinstr(site, "_", "",.)  if inlist(site, "lnu_a", "lnu_b")

save temp, replace

**********************************************
* deal with sites with multiple variants; need to consider these seperately as they have multiple coverage values per site and want to keep all of them for analysis

noi di "where reported as multiple variants"

use temp, clear
keep if strpos(site, "_1") | strpos(site, "_2") | strpos(site, "_3") | strpos(site, "_4") | strpos(site, "_6") 
split site, parse(_)
drop site
destring site2, replace
replace site2 = 10 + site2 if site1=="grla"
replace site1 = "grla/gyra" if inlist(site1, "grla", "gyra")
rename site2 sitenum
rename site1 site
reshape wide value, i(sample method site) j(sitenum)
tab site

* merge with the list of discrepants - this matches all the multiple varient sites
merge 1:1 sample site using discrepancies, update
drop if _merge==1

preserve 
keep if _merge==3
save dis1, replace
restore

* create a list of the sites that haven't been matches to something in multiple variants set
drop if _merge==3
keep sample site  valuegenefinder valuetypewriter valuezam set type ValueAll agree
save dis_remain, replace

***************************************************
* match the single variant sites
noi di "single variants"
use temp, clear
drop if strpos(site, "_1") | strpos(site, "_2") | strpos(site, "_3") | strpos(site, "_4") | strpos(site, "_6") 
* using the leftover from the multiple variants to make sure no site isn't matched
merge 1:1 sample site using dis_remain, update
drop if _merge==1
assert _merge==3
drop _merge

save dis2, replace

***************************************************************************
* Typewriter coverage: Analysis: 
***********************************************************************
*start with single sites as easier
use dis2, clear

* create borderline value flag: <75 = clear no, 75-85 = borderline, 85+ = clear yes 
gen borderline = "no" 
replace borderline = "borderline" if value> 75 | (site=="blaz" & value> 25)
replace borderline= "yes" if value >85 | (site=="blaz" & value > 35)

* here yes = clearly positive, borderline = unsure, no = clearly negative
gen bcount = (borderline=="borderline")
gen ycount = (borderline=="yes")
gen ncount = (borderline=="no")

* creates summary statistics; look at log file for results
summ bcount
noi di r(sum) " out of " r(N) " all sites are borderline +/- 5 from threshold"
noi di r(sum)/r(N) *100 " percent"

foreach value in "AAP" "APA" "APP" "PPA" {
noi di "`value'"

summ bcount if ValueAll=="`value'"
noi di r(sum) " out of " r(N) " sites are borderline +/- 5 from threshold"
noi di r(sum)/r(N) *100 " percent"

summ ycount if ValueAll=="`value'"
noi di r(sum) " out of " r(N) " sites are 5 or more above threshold"
noi di r(sum)/r(N) *100 " percent"

summ ncount if ValueAll=="`value'"
noi di r(sum) " out of " r(N) " sites are 5 or more below threshold"
noi di r(sum)/r(N) *100 " percent"

noi tab site borderline if ValueAll=="`value'"
noi tab set borderline if ValueAll=="`value'"
}


****************************************************************************
* sites with multiple responses

use dis1, clear

* there are two options here: either consider whether any value in the set is borderline, or look at whether the maximum value is borderline.
* given that the variants are later combined to give the overall site prediciton, I'm going to use the maximum
* this creates inaccuracy in chromosonal predictions, as they are based on combonations of genes. However they make up a very small percentage of the errors, so
* I've also eyeballed the data there and there is no clear pattern of borderline results

* find the max across the site variants
gen maxvalue = 0
foreach x in 1 2 3 4 6 11 12 13{
	replace maxvalue= value`x' if value`x'> maxvalue & value`x'!=.
}

* same analysis now as with single sites; see above 
gen borderline = "no" 
replace borderline = "borderline" if maxvalue> 75 | (site=="blaz" & maxvalue> 25)
replace borderline= "yes" if maxvalue >85 | (site=="blaz" & maxvalue > 35)

gen bcount = (borderline=="borderline")
gen ycount = (borderline=="yes")
gen ncount = (borderline=="no")

* summary results; see log file
summ bcount
noi di r(sum) " out of " r(N) " sites (maximised over multiple varients of each site) are borderline +/- 5 from threshold"
noi di r(sum)/r(N) *100 " percent"

foreach value in "AAP" "APA" "APP" "PAP" "PPA" {
noi di "`value'"

summ bcount if ValueAll=="`value'"
noi di r(sum) " out of " r(N) " sites (maximised over multiple varients of each site) are borderline +/- 5 from threshold"
noi di r(sum)/r(N) *100 " percent"

summ ycount if ValueAll=="`value'"
noi di r(sum) " out of " r(N) " sites (maximised over multiple varients of each site) are 5 or more above threshold"
noi di r(sum)/r(N) *100 " percent"

summ ncount if ValueAll=="`value'"
noi di r(sum) " out of " r(N) " sites (maximised over multiple varients of each site) are 5 or more below threshold"
noi di r(sum)/r(N) *100 " percent"

noi tab site borderline if ValueAll=="`value'"
noi tab set borderline if ValueAll=="`value'"
}

* no evidence from Typewriter

exit

* work below was stopped due to the lack of evidence from typewriter making it unlikely much more would be seen in mykrobe results.

*******************************************************************
* MYKROBE: data clean + merge
********************************************************************
noi di _n(5) _dup(80) "=" _n " Mykrobe values" _n _dup(80) "="
********************************************************************

use gene_coverage_zam, clear

* need to clean site name values 
**** match up things I think are the same but with slightly different names
gen site = lower(gene)
replace site =  subinstr(site, "(", "",.)
replace site =  subinstr(site, ")", "",.)
replace site =  subinstr(site, "'", "",.)

replace site = "ant9ib" if site=="ant9_ibspc"
replace site = "ant9ia" if site=="ant9_iaspc"
replace site= "aac6aph2" if strpos(site, "aac")

replace site = subinstr(site, "acr", "arc",.) if strpos(site, "acr")
replace site = "ccrca"  if site == "ccrc_a"
replace site = "ccrcb"  if site == "ccrc_b"
replace site = "ccrcc"  if site == "ccrc_c"
replace site =  subinstr(site, "_", "",.)  if inlist(site, "dfr_a", "dfr_c", "dfr_d", "dfr_g", "dfr_k")
replace site =  subinstr(site, "_", "",.)  if inlist(site, "fus_b", "fus_c")
replace site =  subinstr(site, "_", "",.)  if strpos(site, "qac")

replace site = "grla/gyra" if inlist(site, "grla", "gyra")

duplicates drop sample site percent_coverage median_depth copy_number, force
* Note: there are some duplicated results with different coverage - queried with P but no reponse

* merge in with discrenpancies list; mykrobe doesn't report site varients, so can do single analysis here

merge m:1 sample site using discrepancies, update
drop if _merge==1
assert _merge==3

bysort sample site: gen count=_N
tab count



* create borderline value flag: <75 = clear no, 75-85 = borderline, 85+ = clear yes 
gen borderline = "no" 
replace borderline = "borderline" if percent> 75 | (site=="blaz" & percent> 25)
replace borderline= "yes" if percent >85 | (site=="blaz" & percent > 35)



**************************


