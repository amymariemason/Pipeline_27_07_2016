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
noi di  "Recall order is  Mykrobe Mykrobe Genefinder"


noi tab site ValueA, m
gen agree  = (ValueA=="AAA"| ValueA=="PPP" )
summ agree
noi di r(sum) " out of " r(N) " site predictions agree"
noi di r(sum)/_N*100

noi di "largest disagreements"
noi tab site if agree!=1, sort
noi tab site ValueA if agree!=1

summ agree if strpos(type, "Aqu") |strpos(type, "Chro")
noi di r(sum) " out of " r(N) " aquired or chromosonal resistance site predictions agree"
noi di r(sum)/r(N)*100

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


*sample discrenpancies

noi list sample site Value if agree!=1


******************************************************
noi di _n(5) _dup(80) "=" _n "Site by site predictions" _n _dup(80) "="
**************************

noi tab site ValueAll, m 


noi di _n(5) _dup(80) "=" _n "McNemar relative differences, site by site" _n _dup(80) "="
****************************************
noi di "three ways tables"

use pipeline_clean_all_values_wide, clear
*keep if strpos(type, "viru")
noi bysort site: table valuegenefinder valuetypewriter valuezam 
noi di "site == all"
noi table valuegenefinder valuetypewriter valuezam 

noi di "typewriter vs genefinder"

use pipeline_clean_all_values_wide, clear
*keep if strpos(type, "viru")
table valuegenefinder valuetypewriter 

gen match = ( valuetypewriter== valuegen)

gen AA=( valuegenefinder=="A" &  valuetypewriter=="A")
gen AP=( valuegenefinder=="A" &  valuetypewriter=="P")
gen PA=( valuegenefinder=="P" &  valuetypewriter=="A")
gen PP=( valuegenefinder=="P" &  valuetypewriter=="P")

sort site

* collapse to produce graph
collapse (count) match (sum) AA AP PA PP summatch= match, by(site type)

for any lci uci diff:gen X=.
local max=_N

* create difference intervals

forvalues i=1(1)`max'{
display site[`i']
local a= AA[`i']
local b= AP[`i']
local c= PA[`i']
local d= PP[`i']
mcci `a' `b' `c' `d'
replace diff = r(D_f)  if _n==`i'
replace lci = r(lb_D_f)  if _n==`i'
replace uci = r(ub_D_f)  if _n==`i'
}

* create size of difference ordering for graph
gen prop = abs(diff)
gsort -prop
gen num =_n
gen method= "tg"
labmask num, values(site)

save tgdata, replace


gsort type -prop
gen num2=_n
labmask num2, values(site)

* draw graph
#delimit ;
twoway rcapsym lci uci num2 if method=="tg", s(i) ||
scatter diff num2 if method=="tg" & strpos(type, "Aqu") ,mcolor(red) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0))  ||
scatter diff num2 if method=="tg" & strpos(type, "Chromo") ,mcolor(blue) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) 	||
scatter diff num2 if method=="tg" & strpos(type, "viru") ,mcolor(green) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) ||
scatter diff num2 if method=="tg" & strpos(type, "ccr") ,mcolor(black) msize(small) 
	xlabel(1(1)85, valuelabel angle(90) labsize(tiny))  ylabel(-0.05(0.05)0.05,angle(0)) 	
title("Genefinder vs. Typewriter agreement",size(small) )  
subtitle("More cases where Genefinder=A and Typewriter=P is positive", size(small))
graphregion(fcolor(white)) xtitle("")  
legend(order( 2 3 4 5) label(2 "Aquired Resistance sites") label(3 "Chromosonal mutation sites")
	label(4 "Virulence sites") label(5 "Ccr sites" ));
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\method_disagreements_tg.tif", as(tif) replace


noi di "typewriter vs mykrobe"

use pipeline_clean_all_values_wide, clear
*keep if strpos(type, "viru")
table valuezam valuetypewriter 

gen match = ( valuetypewriter== valuezam)

gen AA=( valuezam=="A" &  valuetypewriter=="A")
gen AP=( valuezam=="A" &  valuetypewriter=="P")
gen PA=( valuezam=="P" &  valuetypewriter=="A")
gen PP=( valuezam=="P" &  valuetypewriter=="P")

sort site

* collapse to produce graph
collapse (count) match (sum) AA AP PA PP summatch= match, by(site type)

for any lci uci diff:gen X=.
local max=_N


forvalues i=1(1)`max'{
display site[`i']
local a= AA[`i']
local b= AP[`i']
local c= PA[`i']
local d= PP[`i']
mcci `a' `b' `c' `d'
replace diff = r(D_f)  if _n==`i'
replace lci = r(lb_D_f)  if _n==`i'
replace uci = r(ub_D_f)  if _n==`i'
}

gen prop = abs(diff)
gsort -prop
gen num =_n
gen method= "tz"
labmask num, values(site)

save tzdata, replace


gsort type - prop
gen num2=_n
labmask num2, values(site)

#delimit ;
twoway rcapsym lci uci num2 if method=="tz", s(i) ||
scatter diff num2 if method=="tz" & strpos(type, "Aqu") ,mcolor(red) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0))  ||
scatter diff num2 if method=="tz" & strpos(type, "Chromo") ,mcolor(blue) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) 	||
scatter diff num2 if method=="tz" & strpos(type, "viru") ,mcolor(green) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) ||
scatter diff num2 if method=="tz" & strpos(type, "ccr") ,mcolor(black) msize(small) 
	xlabel(1(1)85, valuelabel angle(90) labsize(tiny))  ylabel(-0.05(0.05)0.05,angle(0)) 	
title("Mykrobe vs. Typewriter agreement",size(small) )  
subtitle("More cases where Mykrobe=A and Typewriter=P is positive", size(small) )
graphregion(fcolor(white)) xtitle("")  
legend(order( 2 3 4 5) label(2 "Aquired Resistance sites") label(3 "Chromosonal mutation sites")
	label(4 "Virulence sites") label(5 "Ccr sites" ));
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\method_disagreements_tz.tif", as(tif) replace


noi di "mykrobe vs genefinder"

use pipeline_clean_all_values_wide, clear
*keep if strpos(type, "viru")
table valuegenefinder valuezam

gen match = ( valuezam== valuegen)

gen AA=( valuegenefinder=="A" &  valuezam=="A")
gen AP=( valuegenefinder=="A" &  valuezam=="P")
gen PA=( valuegenefinder=="P" &  valuezam=="A")
gen PP=( valuegenefinder=="P" &  valuezam=="P")

sort site

* collapse to produce graph
collapse (count) match (sum) AA AP PA PP summatch= match, by(site type)

for any lci uci diff:gen X=.
local max=_N


forvalues i=1(1)`max'{
display site[`i']
local a= AA[`i']
local b= AP[`i']
local c= PA[`i']
local d= PP[`i']
mcci `a' `b' `c' `d'
replace diff = r(D_f)  if _n==`i'
replace lci = r(lb_D_f)  if _n==`i'
replace uci = r(ub_D_f)  if _n==`i'
}

gen prop = abs(diff)
gsort -prop
gen num =_n
gen method= "zg"
labmask num, values(site)

save zgdata, replace

gsort type -prop
gen num2=_n
labmask num2, values(site)

#delimit ;
twoway rcapsym lci uci num2 if method=="zg", s(i) ||
scatter diff num2 if method=="zg" & strpos(type, "Aqu") ,mcolor(red) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0))  ||
scatter diff num2 if method=="zg" & strpos(type, "Chromo") ,mcolor(blue) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) 	||
scatter diff num2 if method=="zg" & strpos(type, "viru") ,mcolor(green) msize(small) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) ||
scatter diff num2 if method=="zg" & strpos(type, "ccr") ,mcolor(black) msize(small) 
	xlabel(1(1)85, valuelabel angle(90) labsize(tiny))  ylabel(-0.05(0.05)0.05,angle(0)) 	
title("Genefinder vs. Mykrobe agreement",size(small) )  
subtitle("More cases where Genefinder=A and Mykrobe=P is positive", size(small) )
graphregion(fcolor(white)) xtitle("")
legend(order( 2 3 4 5) label(2 "Aquired Resistance sites") label(3 "Chromosonal mutation sites")
	label(4 "Virulence sites") label(5 "Ccr sites" ));
#delimit cr
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\method_disagreements_zg.tif", as(tif) replace



*******
*combine the graphs
*******
noi di "combine graphs"
use tgdata, clear
append using tzdata
append using zgdata


* order sites by abs max diff
sort site
by site: egen ordervalue = total(abs(diff) )
gsort -ordervalue

* create variable to order on graph gsort method -ordervalue2 site
gsort  method type -ordervalue site 
by method: gen newid=_n
order newid
gsort -ordervalue site

labmask newid, values(site)

* create graphs for combo

#delimit ;
twoway rcapsym lci uci newid if method=="zg", s(i) ||
scatter diff newid if method=="zg" & strpos(type, "Aqu") ,mcolor(red) msize(vsmall) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0))  ||
scatter diff newid if method=="zg" & strpos(type, "Chromo") ,mcolor(blue) msize(vsmall) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) 	||
scatter diff newid if method=="zg" & strpos(type, "viru") ,mcolor(green) msize(vsmall) 
	xlabel(1(1)85, valuelabel angle(90))  ylabel(-0.05(0.05)0.05,angle(0)) ||
scatter diff newid if method=="zg" & strpos(type, "ccr") ,mcolor(black) msize(vsmall) 
	xlabel(1(1)85, valuelabel angle(90) labsize(tiny))  ylabel(-0.05(0.05)0.05,angle(0)) 	
title("Genefinder vs. Mykrobe agreement",size(small) )  
subtitle("More cases where Genefinder=P and Mykrobe=A is positive", size(small) )
graphregion(fcolor(white)) xtitle("")
legend(order( 2 3 4 5) label(2 "Aquired Resistance sites") label(3 "Chromosonal mutation sites")
	label(4 "Virulence sites") label(5 "Ccr sites" ))
	fysize(55);
#delimit cr	
graph save graph3, replace


#delimit ;
twoway rcapsym lci uci newid if method=="tg", s(i) ||
scatter diff newid if method=="tg" & strpos(type, "Aqu") ,mcolor(red) msize(vsmall) 
	xlabel(none) xtitle("")  ylabel(-0.05(0.05)0.05,angle(0))  ||
scatter diff newid if method=="tg" & strpos(type, "Chromo") ,mcolor(blue) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-0.05(0.05)0.05,angle(0)) 	||
scatter diff newid if method=="tg" & strpos(type, "viru") ,mcolor(green) msize(vsmall) 
	xlabel(none) xtitle("")  ylabel(-0.05(0.05)0.05,angle(0)) ||
scatter diff newid if method=="tg" & strpos(type, "ccr") ,mcolor(black) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-0.05(0.05)0.05,angle(0)) 	
title("Genefinder vs. Typewriter agreement",size(small) )  
subtitle("More cases where Genefinder=P and Typewriter=A is positive", size(small) )
graphregion(fcolor(white)) xtitle("") legend(off) fysize(30);
#delimit cr
graph save graph2, replace


#delimit ;
twoway rcapsym lci uci newid if method=="tz", s(i) ||
scatter diff newid if method=="tz" & strpos(type, "Aqu") ,mcolor(red) msize(vsmall) 
	xlabel(none) xtitle("")  ylabel(-0.05(0.05)0.05,angle(0))  ||
scatter diff newid if method=="tz" & strpos(type, "Chromo") ,mcolor(blue) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-0.05(0.05)0.05,angle(0)) 	||
scatter diff newid if method=="tz" & strpos(type, "viru") ,mcolor(green) msize(vsmall) 
	xlabel(none) xtitle("")  ylabel(-0.05(0.05)0.05,angle(0)) ||
scatter diff newid if method=="tz" & strpos(type, "ccr") ,mcolor(black) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-0.05(0.05)0.05,angle(0)) 	
title("Mykrobe vs. Typewriter agreement",size(small) )  
subtitle("More cases where Mykrobe=P and Typewriter=A is positive", size(small) )
graphregion(fcolor(white)) xtitle("")  
legend(off)
fysize(30);
#delimit cr
graph save graph1, replace


graph combine graph1.gph graph2.gph graph3.gph, cols(1) ycommon

graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\method_disagreements_all.tif", as(tif) width(2550) replace

graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\method_disagreements_all.eps", as(eps) replace



