* Look at the coverage of the sites that were discrepancies
* uses from typewriter and mykrobe to do this

set li 130

cap log close
log using site.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

use pipeline_clean_all_values_wide, clear

**************************
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
* hypothesis: this is due to borderline coverage

use coverage_tw, clear


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
* deal with sites with multiple variants
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


merge 1:1 sample site using discrepancies, update
drop if _merge==1

preserve 
keep if _merge==3
save dis1, replace
restore

drop if _merge==3
keep sample site  valuegenefinder valuetypewriter valuezam set type ValueAll agree
save dis_remain, replace

***************************************************
* deal with rest
noi di "single variants"
use temp, clear
drop if strpos(site, "_1") | strpos(site, "_2") | strpos(site, "_3") | strpos(site, "_4") | strpos(site, "_6") 
merge 1:1 sample site using dis_remain, update
drop if _merge==1
assert _merge==3
drop _merge

save dis2, replace

***************************************************************************
* Analysis: start with single sites as easier

use dis2, replace

* create borderline value flag: <75 = clear no, 75-85 = borderline, 85+ = clear yes 
gen borderline = "no" 
replace borderline = "borderline" if value> 75 | (site=="blaz" & value> 25)
replace borderline= "yes" if value >85 | (site=="blaz" & value > 35)

gen bcount = (borderline=="borderline")
gen ycount = (borderline=="yes")
gen ncount = (borderline=="no")

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

* Not clear to me that borderline values are the problem. 
*******************

*