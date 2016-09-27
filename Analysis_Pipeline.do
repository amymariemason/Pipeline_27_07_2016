* Compare the three methods site by site

set li 130

cap log close
log using analysis1.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

*************************************************
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets
use pipeline_clean_all_values_wide, clear

****************************************

******************************************************
noi di _n(5) _dup(80) "=" _n " Overall summary" _n _dup(80) "="
**************************
noi di  "Recall order is  Typewriter Mykrobe Genefinder"
noi tab ValueAll, sort
gen agree =inlist(ValueAll, "AAA", "PPP")
summ agree 
noi di r(sum) " out of " r(N) " predictions are identical between all three methods"

noi tab ValueAll if strpos(type, "Aquired"), sort
summ agree  if strpos(type, "Aquired")
noi di r(sum) " out of " r(N) " predictions are identical between all three methods on aquired resistance prediction"

noi tab ValueAll if strpos(type, "Chrom"), sort
summ agree  if strpos(type, "Chrom")
noi di r(sum) " out of " r(N) " predictions are identical between all three methods on chromosonal resistance prediction"

noi tab ValueAll if strpos(type, "ccr"), sort
summ agree  if strpos(type, "ccr")
noi di r(sum) " out of " r(N) " predictions are identical between all three methods on ccr prediction"

noi tab ValueAll if strpos(type, "viru"), sort
summ agree  if strpos(type, "viru")
noi di r(sum) " out of " r(N) " predictions are identical between all three methods on virulence prediction"


noi tab ValueAll if strpos(type, "Chrom")|strpos(type, "Aquired"), sort
summ agree  if strpos(type, "Chrom")|strpos(type, "Aquired")
noi di r(sum) " out of " r(N) " predictions are identical between all three methods on antibiotic prediction"


******************************************************
noi di _n(5) _dup(80) "=" _n "Site by site predictions" _n _dup(80) "="
**************************

noi tab site ValueAll, m 
