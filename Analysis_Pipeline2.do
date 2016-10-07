* Compare the three methods in phenotype predictions

set li 130

cap log close
log using analysis2.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

*********** create phenotype predictions by site
use anti_panel_all, clear

noi di" make bar chart of which methods differ on which antibiotic"
preserve
contract  site  value*
rename valueall all
drop value*

reshape wide _freq, i(site) j(all) string
rename _freq* *


* bar graph of predictions (antibiotics)
graph bar (asis)  rrr rrs rsr rss srs sss,  over(site, label(angle(90))) stack title("Combinations of Results") subtitle("Results given as Genefinder Mykrobe Typewriter ") legend( label(1 "rrr") label( 2 "rrs" ) label (3 "rsr") label (4 "rss" ) label (5 "srs") label( 6 "sss" )) bar(1, color(gs8)) bar(2, color(cyan)) bar(3, color(red)) bar(4, color(blue)) bar(5, color(yellow)) bar(6, color(gs12))
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\phenotype_disagreements.tif", as(tif) width(2550) replace

restore
*************************************************************

noi di "compare to gold standard"

* in table
noi tab gold valueall

* in graph
preserve
contract  site  gold valueall

reshape wide _freq, i(site gold) j(valueall) string
rename _freq* *


gen graphx = site + " " + gold
drop if !inlist(gold, "R", "S")
sort gold site
 gen num =_n
 replace num=num+1 if gold=="S"
labmask num, values(site)
gen goldlabel = "Gold standard " + gold

* bar graph of predictions (antibiotics)
graph bar (asis) rrr rrs rsr rss srs sss if inlist(golds,"R", "S"),graphregion(color(white))    over(site, label(angle(90) labsize(tiny)) )  over(goldlabel, label(labsize(small)) ) stack  title("Prediction combinations") subtitle("Results given as Genefinder Mykrobe Typewriter") legend(rows(2) label(1 "rrr") label( 2 "rrs" ) label (3 "rsr") label (4 "rss" ) label (5 "srs") label( 6 "sss" )) bar(1, color(gs8)) bar(2, color(cyan)) bar(3, color(red)) bar(4, color(blue)) bar(5, color(yellow)) bar(6, color(gs12)) ylabel(0(500)1000 1400)
graph export "E:\users\amy.mason\Pipeline_27_07_2016\Graphs_Outputs\phenotype_disagreements_withgold.tif", as(tif) width(2550) replace

restore



* sensitivity/ specificity







