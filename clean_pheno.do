* clean pipeline phenotype data

set li 130

cap log close
log using clean2.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets


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
 

************************************************
noi di _n(5) _dup(80) "=" _n " 2 Use PCR/Array to interpret phenotypes" _n _dup(80) "="

************************************************

noi dis "based on entero pcr/array"
foreach x in A B C D E G H I J LR P U{
	gen SE`x' = ""
	replace SE`x' = "R" if strpos( entero_pcr, "`x'") | strpos(entero_array, "`x'")
	replace SE`x' = "S" if !strpos( entero_pcr, "`x'") & !strpos(entero_array, "`x'") & (entero_pcr!=""|entero_array!="")
	replace SE`x' = "NT" if strpos( lower(entero_pcr), "nt") & strpos(lower(entero_array), "nt")
	noi tab SE`x', m		
			}
		
noi di "based on exfol_pcr/ array"
		foreach x in A B D {
	gen ET`x' = ""
	replace ET`x' = "R" if strpos( exfol_pcr, "`x'") | strpos(exfol_array, "`x'")
	replace ET`x' = "S" if !strpos( exfol_pcr, "`x'") & !strpos(exfol_array, "`x'")  & (exfol_pcr!=""|exfol_array!="")
	replace ET`x' = "NT" if strpos( lower(exfol_pcr), "nt") & strpos(lower(exfol_array), "nt")
	noi tab ET`x', m		
			}
				
noi di "based on tst pcr/array"
		
	gen TSST1 = "S" if strpos( tst_pcr, "-") | strpos(tst_array, "-")
	replace TSST1 = "R" if strpos( tst_pcr, "+") | strpos(tst_array, "+")
	replace TSST1 = "NT" if strpos( lower(tst_pcr), "nt") & strpos(lower(tst_array), "nt")
noi tab TSST1, m

noi di "based on PVL pcr and array"	
	gen PVL = "S" if strpos( PVLbyPCR, "-") | strpos(PVL_array, "-")
	replace PVL = "R" if strpos( PVLbyPCR, "+") | strpos(PVL_array, "+")
	replace PVL = "NT" if strpos( lower(PVLbyPCR), "nt") & strpos(lower(PVL_array), "nt")
noi tab PVL,m 

noi di "based on MecA pcr"	
	gen MECA = "S" if  mecA_pcr=="-"
	replace MECA = "R" if mecA_pcr=="+"
	replace MECA = "NT" if mecA_pcr=="NT"
noi tab MECA, m	

noi di "based on MecC pcr"	
	gen MECC = "S" if  mecC_pcr=="-"
	replace MECC = "R" if mecC_pcr=="+"
	replace MECC = "NT" if mecC_pcr=="NT"
noi tab MECC, m
		
drop  comments mecA_pcr mecC_pcr PVLbyPCR PVL_array entero_pcr entero_array exfol_pcr exfol_array tst_pcr tst_array
compress

* reduce variation in answers 
************************************************
noi di _n(5) _dup(80) "=" _n "2 reduce variation in notation from method to method" _n _dup(80) "="

************************************************


foreach k in Penicillin Methicillin Ciprofloxacin Erythromycin Clindamycin Tetracycline Fusidicacid Gentamicin Rifampicin Trimethoprim Vancomycin Mupirocin ClindDtest SEA SEB SEC SED SEE SEG SEH SEI SEJ SELR SEP SEU ETA ETB ETD TSST1 PVL MECA MECC{
		noi di "`k'"
		gen ambi = 1 if inlist(upper(`k'), "ND", "NT", "LOW COV", "MIXED")
		replace ambi=1 if inlist(upper(`k'), "NOT DONE", "N/A", "NF", "NA")
		summ ambi
		noi di r(sum) " values are in list ND, NT, LOW COV, MIXED, N/A, NOT DONE, replacing with NA"
		replace `k' = "NA" if ambi==1
		drop ambi
		
		gen ambi = 1 if `k'=="B"
		summ ambi
		noi di r(sum)  " reported as B; replace with I for consistancy between methods "
		replace `k' ="I" if ambi==1
		drop ambi
		
		gen ambi = 1 if `k'==""
		summ ambi
		noi di r(sum)  " entirely missing values; replace with -"
		replace `k' ="-" if ambi==1
		drop ambi
		
		noi di "final version"
		noi tab `k', m
		}

		
		

	
* fix clindamycin tests
************************************************
noi di _n(5) _dup(80) "=" _n "3   use D-test to fix clindamycin test results" _n _dup(80) "="

************************************************

gen count =1 if  Clindamycin=="S" &  ClindDtest == "POS" 
summ  count if count==1
noi di r(N) " clinddtest results are positive when clindamycin results = S -> replace clindamycin test with R"
replace  Clindamycin= "R" if count == 1
noi tab Clindamycin, m
drop *ClindDtest count

	noi di "SAVE: one record per sample with gold standard profile"
	save pipeline_gold_clean_wide, replace

		
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
	
	
noi cd E:\users\amy.mason\Pipeline_27_07_2016\	