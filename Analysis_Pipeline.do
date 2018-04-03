************************************************
* ANALYSIS_PIPELINE.DO 
************************************************

* Analyses the three methods on a site by site basis; creates summary tables agraphs

*Inputs:  pipeline_clean_all_values_wide (from clean_predict.do)
* Outputs :  method_disagreements_all.wmf ( graph comparing the differences in the methods in three pairs: typewriter vs. genefinder, genefinder vs mykrobe, and mykrobe vs typewriter
* graphs saved to Graph folder
* Written by: Amy Mason


set li 130

cap log close
log using analysis1.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Datasets"

*************************************************
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Datasets"
use pipeline_clean_all_values_wide, clear

****************************************

******************************************************
noi di _n(5) _dup(80) "=" _n " Overall summary" _n _dup(80) "="
******************************************************
* tables summarizing the % agreements between methods
* look at log file to read summary

noi di  "Recall order is  Genefinder Mykrobe Typewriter"


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


*which samples have the most discrenpacies
noi di "site discrenpacies by sample"
noi tab sample if agree!=1, sort
noi list sample site Value if agree!=1


*which set  have the most discrenpacies
noi di "site discrenpacies by set"
noi tab set if agree!=1, sort
noi tab set agree, row
noi tab set agree, chi2
noi list set site Value if agree!=1

*which set and type  have the most discrenpacies
noi di "site discrenpacies by set and type"
gen new = set + type
noi tab new if agree!=1, sort
noi tab new agree, row
noi tab new agree, chi2

*Kappa statistics for paper
noi di "Kappa statistics site by site"
gen gf = ( valuegenefinder=="P")
gen tw = ( valuet=="P")
gen mk = ( valuez=="P")

noi di "genefinder compared with typewriter"
noi kap  gf tw
noi di "typewriter compared to mykrobe"
noi kap tw mk
noi di "mykrobe compared to genefinder"
noi kap mk gf

noi di "all three"
noi kap tw mk gf

drop gf tw mk


*********************************************************
noi di _n(5) _dup(80) "=" _n "Site by site predictions" _n _dup(80) "="
*********************************************************
* which sites have the most discrepancies, also table for paper
noi tab site ValueAll, m 


noi di _n(5) _dup(80) "=" _n "McNemar relative differences, site by site" _n _dup(80) "="
********************************************************************************
noi di "three ways tables"

use pipeline_clean_all_values_wide, clear
*keep if strpos(type, "viru")
noi bysort site: table valuegenefinder valuetypewriter valuezam 
noi di "site == all"
noi table valuegenefinder valuetypewriter valuezam 

***********************************************************************
* Compare typewriter with genefinder
**********************************************************************


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

* calculate mcnemar chi, is there a significant bias in which was discrepancies go
preserve

noi di "McNemar chi"

collapse (count) match (sum) AA AP PA PP summatch= match

for any McNchi pvalue chi95:gen X=.

replace McNchi = ((abs(AP - PA)-1)^2)/(AP+PA)
replace pvalue = 1- chi2(1,McNchi)

replace chi95 = 3.841
gen sig = "Significant" if McNchi> chi95 & McNchi!=.

noi list McNchi pvalue sig if sig!=""

restore

gen GY = ( valuegenefinder=="P")
gen TY = ( valuet=="P")
gen ZY = ( valuez=="P")

mcc GY TY

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

* save graph
graph export "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Graphs_Outputs\method_disagreements_tg.wmf", as(wmf) replace

exit
***********************************************************************
* Compare mykrobe with typewriter
**********************************************************************

noi di "typewriter vs mykrobe"

use pipeline_clean_all_values_wide, clear
* summary table of results
table valuezam valuetypewriter 
* totals for calculating mcnemar
gen match = ( valuetypewriter== valuezam)
gen AA=( valuezam=="A" &  valuetypewriter=="A")
gen AP=( valuezam=="A" &  valuetypewriter=="P")
gen PA=( valuezam=="P" &  valuetypewriter=="A")
gen PP=( valuezam=="P" &  valuetypewriter=="P")
sort site

* calculate mcnemar chi, is there a significant bias in results
preserve
noi di "McNemar chi"
collapse (count) match (sum) AA AP PA PP summatch= match
for any McNchi pvalue chi95:gen X=.
replace McNchi = ((abs(AP - PA)-1)^2)/(AP+PA)
replace pvalue = 1- chi2(1,McNchi)
replace chi95 = 3.841
gen sig = "Significant" if McNchi> chi95 & McNchi!=.
noi list McNchi pvalue sig if sig!=""
restore

* collapse to produce graph
collapse (count) match (sum) AA AP PA PP summatch= match, by(site type)

* create difference confidence intervals
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

* order by size of difference
gen prop = abs(diff)
gsort -prop
gen num =_n
gen method= "tz"
labmask num, values(site)

save tzdata, replace

* make ordered labels for graph; inverse size of difference
gsort type - prop
gen num2=_n
labmask num2, values(site)
* make graph
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
* save graph
graph export "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Graphs_Outputs\method_disagreements_tz.wmf", as(wmf) replace

***********************************************************************
* Compare mykrobe with genefinder
**********************************************************************
noi di "mykrobe vs genefinder"


use pipeline_clean_all_values_wide, clear
* create summary table
table valuegenefinder valuezam
* create markers to total for mcnemar
gen match = ( valuezam== valuegen)

gen AA=( valuegenefinder=="A" &  valuezam=="A")
gen AP=( valuegenefinder=="A" &  valuezam=="P")
gen PA=( valuegenefinder=="P" &  valuezam=="A")
gen PP=( valuegenefinder=="P" &  valuezam=="P")

sort site

* calculate mcnemar chi, is there a significant bias in results
preserve

noi di "McNemar chi"

collapse (count) match (sum) AA AP PA PP summatch= match

for any McNchi pvalue chi95:gen X=.

replace McNchi = ((abs(AP - PA)-1)^2)/(AP+PA)
replace pvalue = 1- chi2(1,McNchi)

replace chi95 = 3.841
gen sig = "Significant" if McNchi> chi95 & McNchi!=.

noi list McNchi pvalue sig if sig!=""



restore

* collapse instead  to produce site by site graph
collapse (count) match (sum) AA AP PA PP summatch= match, by(site type)

for any lci uci diff:gen X=.
local max=_N

* calculate confidence intervals for difference
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

* order in by size of difference
gen prop = abs(diff)
gsort -prop
gen num =_n
gen method= "zg"
labmask num, values(site)

save zgdata, replace
* order by reverse difference
gsort type -prop
gen num2=_n
labmask num2, values(site)

* make graph

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
*save graph
graph export "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Graphs_Outputs\method_disagreements_zg.wmf", as(wmf) replace



********************************************************************
*combine the three comparision graphs
********************************************************************
noi di "combine graphs"
* import graph data from above
use tgdata, clear
append using tzdata
append using zgdata


* order sites by abs max diff
sort site
by site: egen ordervalue = total(abs(diff) )
gsort -ordervalue

* create variable to order on graph 
gsort  method type -ordervalue site 
by method: gen newid=_n
order newid
gsort -ordervalue site

labmask newid, values(site)

* change values to percentage (as Sarah's request)

replace diff= diff*100
replace lci = lci *100
replace uci = uci *100

* create graphs for combining so that orders match

#delimit ;
twoway rcapsym lci uci newid if method=="zg", s(i) ||
scatter diff newid if method=="zg" & strpos(type, "Aqu") ,mcolor(red) msize(vsmall) 
	xlabel(1(1)83, valuelabel angle(90))  ylabel(-5(5)5,angle(0))  ||
scatter diff newid if method=="zg" & strpos(type, "Chromo") ,mcolor(blue) msize(vsmall) 
	xlabel(1(1)83, valuelabel angle(90))  ylabel(-5(5)5,angle(0)) 	||
scatter diff newid if method=="zg" & strpos(type, "viru") ,mcolor(green) msize(vsmall) 
	xlabel(1(1)83, valuelabel angle(90))  ylabel(-5(5)5,angle(0)) ||
scatter diff newid if method=="zg" & strpos(type, "ccr") ,mcolor(black) msize(vsmall) 
	xlabel(1(1)83, valuelabel angle(90) labsize(tiny))  ylabel(-5(5)5,angle(0)) 	
title("Genefinder vs. Mykrobe agreement",size(small) ) 
	graphregion(fcolor(white)) xtitle("Gene")
legend(order( 2 3 4 5) label(2 "Genes with aquired resistance") label(3 "Genes with chromosonal mutation")
	label(4 "Genes with virulence") label(5 "Genes with CCR" ))
	fysize(55)
ytitle("Percentage", size(small))	;
#delimit cr	
*subtitle("More cases where Genefinder=P and Mykrobe=A is positive", size(small) )
graph save graph3, replace


#delimit ;
twoway rcapsym lci uci newid if method=="tg", s(i) ||
scatter diff newid if method=="tg" & strpos(type, "Aqu") ,mcolor(red) msize(vsmall) 
	xlabel(none) xtitle("")  ylabel(-5(5)5,angle(0)) ||
scatter diff newid if method=="tg" & strpos(type, "Chromo") ,mcolor(blue) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-5(5)5,angle(0)) 	||
scatter diff newid if method=="tg" & strpos(type, "viru") ,mcolor(green) msize(vsmall) 
	xlabel(none) xtitle("")  ylabel(-5(5)5,angle(0))||
scatter diff newid if method=="tg" & strpos(type, "ccr") ,mcolor(black) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-5(5)5,angle(0))	
title("Genefinder vs. Typewriter agreement",size(small) )  
ytitle("Percentage", size(small))
graphregion(fcolor(white)) xtitle("") legend(off) fysize(30);
#delimit cr
graph save graph2, replace
*subtitle("More cases where Genefinder=P and Typewriter=A is positive", size(small) )

#delimit ;
twoway rcapsym lci uci newid if method=="tz", s(i) ||
scatter diff newid if method=="tz" & strpos(type, "Aqu") ,mcolor(red) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-5(5)5,angle(0)) ||
scatter diff newid if method=="tz" & strpos(type, "Chromo") ,mcolor(blue) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-5(5)5,angle(0))	||
scatter diff newid if method=="tz" & strpos(type, "viru") ,mcolor(green) msize(vsmall) 
	xlabel(none) xtitle("")  ylabel(-5(5)5,angle(0)) ||
scatter diff newid if method=="tz" & strpos(type, "ccr") ,mcolor(black) msize(vsmall) 
	xlabel(none) xtitle("") ylabel(-5(5)5,angle(0))	
title("Mykrobe vs. Typewriter agreement",size(small) )  
ytitle("Percentage", size(small))
graphregion(fcolor(white)) xtitle("")  
legend(off)
fysize(30);
#delimit cr
graph save graph1, replace
*subtitle("More cases where Mykrobe=P and Typewriter=A is positive", size(small) )
* combine the graphs

graph combine graph1.gph graph2.gph graph3.gph, cols(1) ycommon

*save the graph

graph export "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Graphs_Outputs\method_disagreements_all.wmf", as(wmf) replace

graph export "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Graphs_Outputs\method_disagreements_all.eps", as(eps) replace

graph export "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Graphs_Outputs\method_disagreements_all_0205.tif", as(tif) width(2880) replace

