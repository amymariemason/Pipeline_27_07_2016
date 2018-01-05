* clean pipeline phenotype data

set li 130

cap log close
log using clean2.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd "\\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\Datasets"


*************************************************

use pipeline_data_gs, clear
* to distinguish blank responses in pcrs checked from blacnk boxes formed by long -> wide

******************************************************
noi di _n(5) _dup(80) "=" _n " 1 drop contaminated samples" _n _dup(80) "="
**************************

gen contam=.
replace contam=1 if inlist(sample,"C00012813", "C00001215", "C00001249", "C00012746", "C00012791")
replace contam=1 if inlist(sample,"H113120068-138-1","H121000461-388-2","H131100031-408-2","H111840168-459-2","H111200061-196-2")
bysort sample: gen count=1 if _n==1
summ contam if count==1
noi di r(sum) " samples contaminated"
noi drop if contam==1
drop contam count


******************************************************
noi di _n(5) _dup(80) "=" _n " 2 Match up names from method to method" _n _dup(80) "="
*********************************************************

local list1 "ClindDtest Clindamycin Penicillin Methicillin Ciprofloxacin  Clindamycin  Erythromycin  Tetracycline Fusidicacid Gentamicin Rifampicin Trimethoprim Vancomycin Mupirocin"
local list2 "Dtest Clindamycindisc PEN OXA CIP CLI ERY TET FUS GEN RIF TMP VAN MUP"
local n1 : word count `list1'
local n2 : word count `list2'
if `n1' != `n2' {
        di as err "different counts"
}

forval i = 1/`n1' {
        local first : word `i' of `list1'
        local second : word `i' of `list2'
        noi display " match `first' and `second'"
		assert !(`first' !="" & `second' !="")
		replace `first' = `second' if `first' ==""
		replace `first' = upper(`first')
		drop `second'
		noi tab `first', m
		}
 

 * reduce variation in answers 
************************************************
noi di _n(5) _dup(80) "=" _n "2 reduce variation in notation from method to method (antibiotics) " _n _dup(80) "="

************************************************
local NAcount=0
local Bcount =0
local blankcount=0


foreach k in Penicillin Methicillin Ciprofloxacin Erythromycin Clindamycin Tetracycline Fusidicacid Gentamicin Rifampicin Trimethoprim Vancomycin Mupirocin ClindDtest mecA_pcr mecC_pcr PVLbyPCR PVL_array entero_pcr entero_array exfol_pcr exfol_array tst_pcr tst_array{
		noi di "`k'"
		gen ambi = 1 if inlist(upper(`k'), "ND", "NT", "LOW COV", "MIXED")
		replace ambi=1 if inlist(upper(`k'), "NOT DONE", "N/A", "NF", "NA")
		summ ambi
		noi di r(sum) " values are in list ND, NT, LOW COV, MIXED, N/A, NOT DONE, replacing with NA"
		local NAcount = `NAcount'+r(sum)
		replace `k' = "NA" if ambi==1
		drop ambi
		
		gen ambi = 1 if `k'=="B"
		summ ambi
		noi di r(sum)  " reported as B; replace with I for consistancy between methods "
		local Bcount = `Bcount'+r(sum)
		replace `k' ="I" if ambi==1
		drop ambi
		
		gen ambi = 1 if `k'==""
		summ ambi
		noi di r(sum)  " entirely missing values; replace with NA"
		local blankcount = `blankcount'+r(sum)
		replace `k' ="NA" if ambi==1
		drop ambi
		
		}
noi di "in total " `NAcount' " values are in list ND, NT, LOW COV, MIXED, N/A, NOT DONE, replacing with NA"
noi di "in total " `Bcount' " reported as B; replace with I for consistancy between methods "
noi di "in total " `blankcount' " entirely missing values; replace with NA"
 
 
************************************************
noi di _n(5) _dup(80) "=" _n " 3 Use PCR/Array to interpret phenotypes" _n _dup(80) "="

************************************************

noi dis "based on entero pcr/array"
foreach x in A B C D E G H I J LR P U{
	gen SE`x' = ""
	replace SE`x' = "R" if (strpos( entero_pcr, "`x'") & entero_pcr!="NA" ) | (strpos(entero_array, "`x'")& entero_array!="NA")
	replace SE`x' = "NA" if (strpos( entero_pcr, "NA") & strpos(entero_array, "NA")) 
	replace SE`x' = "S" if !inlist(SE`x', "R", "NA")
	noi tab SE`x', m		
			}
		
noi di "based on exfol_pcr/ array"
		foreach x in A B D {
	gen ET`x' = ""
	replace ET`x' = "R" if (strpos( exfol_pcr, "`x'") & exfol_pcr!="NA") | (strpos(exfol_array, "`x'")& exfol_array!="NA")
	replace ET`x' = "NA" if (strpos( exfol_pcr, "NA") & strpos(exfol_array, "NA"))
	replace ET`x' = "S" if  !inlist(ET`x', "R", "NA")
	noi tab ET`x', m		
			}
				
noi di "based on tst pcr/array"
		
	gen TSST1 = "R" if strpos( tst_pcr, "+") | strpos(tst_array, "+")
	replace TSST1 = "NA" if strpos(tst_pcr, "NA") & strpos(tst_array, "NA") 
	replace TSST1 = "S" if  !inlist(TSST1, "R", "NA")
noi tab TSST1, m

noi di "based on PVL pcr and array"	
	gen PVL = "R" if strpos( PVLbyPCR, "+") | strpos(PVL_array, "+")
	replace PVL = "NA" if ( strpos(PVLbyPCR, "NA") & strpos(PVL_array, "NA"))
	replace PVL = "S" if !inlist(PVL, "R", "NA")
noi tab PVL,m 

noi di "based on MecA pcr"	
	gen MECA = "R" if mecA_pcr=="+"
	replace MECA = "NA" if mecA_pcr=="NA"
	replace MECA = "S" if !inlist(MECA, "R", "NA")
noi tab MECA, m	

noi di "based on MecC pcr"	
	gen MECC = "R" if mecC_pcr=="+"
	replace MECC = "NA" if mecC_pcr=="NA"
	replace MECC = "S" if !inlist(MECC, "R", "NA")
	
noi tab MECC, m
		
drop  comments mecA_pcr mecC_pcr PVLbyPCR PVL_array entero_pcr entero_array exfol_pcr exfol_array tst_pcr tst_array
compress


		

	
* fix clindamycin tests
************************************************
noi di _n(5) _dup(80) "=" _n "4   use D-test to fix clindamycin test results" _n _dup(80) "="

************************************************

gen count =1 if  Clindamycin=="S" &  ClindDtest == "POS" 
summ  count if count==1
noi di r(N) " clinddtest results are positive when clindamycin results = S -> replace clindamycin test with R"
replace  Clindamycin= "R" if count == 1
noi tab Clindamycin, m
drop *ClindDtest count

	noi di "SAVE: one record per sample with gold standard profile"
	save pipeline_gold_clean_wide, replace


**********************************************************************
* save results by site
************************************************
noi di _n(5) _dup(80) "=" _n "5  final test results" _n _dup(80) "="

************************************************

foreach k in Penicillin Methicillin Ciprofloxacin Erythromycin Clindamycin Tetracycline Fusidicacid Gentamicin Rifampicin Trimethoprim Vancomycin Mupirocin SEA SEB SEC SED SEE SEG SEH SEI SEJ SELR SEP SEU ETA ETB ETD TSST1 PVL MECA MECC{
		noi di "`k'"
		noi tab `k', m	
			
}

	
* reshape back and match up names

	rename * value*
	rename valuesample sample
	rename valueset set
	rename valuemethod method

	reshape long value, i(sample) j(site) string
	replace site = lower(site)
	
	

	rename value goldstandard
	
	
	noi di "SAVE: one record per sample-site with gold standard profile"
	save pipeline_gold_clean_long, replace
	
	
noi cd \\me-filer1\home$\am2609\My Documents\MMM work\Programs\Pipeline\	