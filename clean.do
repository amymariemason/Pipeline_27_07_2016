* clean pipeline data

set li 130

cap log close
log using clean.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets
***********************************************************
* 1) Remove samples that are contaminated for all methods
*********************************************************

noi di _n(5) _dup(80) "=" _n " 1 Remove samples that were contaminated for all methods" _n _dup(80) "="

use pipeline_all, clear
preserve
reshape wide value  comment, i(set sample site) j(method) string
noi display _N " sample sites"
restore
gen contam=.
replace contam=1 if inlist(sample,"C00012813", "C00001215", "C00001249", "C00012746", "C00012791")
replace contam=1 if inlist(sample,"H113120068-138-1","H121000461-388-2","H131100031-408-2","H111840168-459-2","H111200061-196-2")
preserve
keep if contam==1
reshape wide value  comment, i(set sample site) j(method) string
* count
bysort sample: gen count=1 if _n==1
summ contam if count==1
noi di r(sum) " samples contaminated"

keep sample set site value* comment* 
compress
sort sample site
order sample site value* comment*

noi save problem_samples, replace


restore
summ contam
noi di r(sum) " values for contaminated samples dropped from study"
noi drop if contam==1
drop contam
compress
save pipeline_noncontam, replace
*****************************************************
* clean site names
******************************************************
use pipeline_noncontam, clear
noi di _n(5) _dup(80) "=" _n "consolidate names/ spelling errors for various sites" _n _dup(80) "=" 

**** match up things I think are the same but with slightly different names
replace site = "ant9ib" if site=="ant9_ibspc"
replace site = "ant9ia" if site=="ant9_iaspc"
replace site= "aac6aph2" if strpos(site, "aac")
replace site= "clinddtest" if strpos(site, "dtest")
replace site="clindamycin" if site=="clindamycindisc"
replace site = subinstr(site, "acr", "arc",.) if strpos(site, "acr")
replace site = "mlst" if site=="st"
replace site = "ccrca"  if site == "ccrc_a"
replace site = "ccrcb"  if site == "ccrc_b"
replace site = "ccrcc"  if site == "ccrc_c"
replace site =  subinstr(site, "_", "",.) if strpos(site, "dfr")
replace site =  subinstr(site, "_", "",.) if strpos(site, "fus")
replace site =  subinstr(site, "_", "",.) if strpos(site, "ile")
replace site =  subinstr(site, "_", "",.) if strpos(site, "lnu")
replace site =  subinstr(site, "_", "",.) if strpos(site, "qac")
replace site =  subinstr(site, "_", "",.) if strpos(site, "van")

***** match up phenotype names

replace site="methicillin" if site=="oxa"
replace site= "ciprofloxacin" if inlist(site, "cip")
replace site= "fusidicacid" if inlist(site, "fus")
replace site= "clindamycin" if inlist(site, "cli")
replace site= "erythromycin" if inlist(site, "ery")
replace site= "gentamicin" if inlist(site, "gen")
replace site= "mupirocin" if inlist(site, "mup")
replace site= "penicillin" if inlist(site, "pen")
replace site= "rifampicin" if inlist(site, "rif")
replace site= "tetracycline " if inlist(site, "tet")
replace site= "trimethoprim" if inlist(site, "tmp")
replace site= "vancomycin" if inlist(site, "van")

replace site = subinstr(site, " ", "",.)
compress
save temp, replace



* some sites split into numbers variants; remerge

noi di "combining numbered variants into single sites"

use temp, clear

drop comments 
keep if strpos(site, "erma") | strpos(site,"ccr") | strpos(site, "se")
drop if inlist(site, "ccrcb", "ccrca", "ccrcc", "seb", "sec", "sed")
reshape wide value, i(set sample method) j(site) string

*erma
replace valueerma="A" if valueerma_1=="A" & valueerma_2=="A" &valueerma==""
replace valueerma="P" if (valueerma_1=="P" | valueerma_2=="P" ) & valueerma==""
assert valueerma!=""
drop valueerma_1 valueerma_2

*ccra
replace valueccra="A" if valueccra=="" & valueccra_1=="A" & valueccra_2=="A" & valueccra_3=="A" & valueccra_4=="A"
replace valueccra="P" if valueccra=="" & (valueccra_1=="P"|  valueccra_2=="P"| valueccra_3=="P"| valueccra_4=="P")
assert valueccra!=""
drop valueccra_*

* ccrb
replace valueccrb="A" if valueccrb=="" & valueccrb_1=="A" & valueccrb_2=="A" & valueccrb_3=="A" & valueccrb_4=="A" & valueccrb_6=="A" 
replace valueccrb="P" if valueccrb=="" & (valueccrb_1=="P" | valueccrb_2=="P"|valueccrb_3=="P"|valueccrb_4=="P"|valueccrb_6=="P" )
* two values blank  - C00011656 & C00011655
assert valueccrb!="" if !strpos(sample, "C0001165")
drop valueccrb_*

* sea
replace valuesea="A" if valuesea_1=="A" & valuesea_2=="A"  & valuesea=="" 
replace valuesea="P" if valuesea_1=="P" | valuesea_2=="P"  & valuesea=="" 
* two values blank  - C00011656 & C00011655
assert valuesea!=""  if !strpos(sample, "C0001165")
drop valuesea_*

*seh
replace valueseh="A" if valueseh_1=="A" & valueseh_2=="A"  & valueseh=="" 
replace valueseh="P" if valueseh_1=="P" | valueseh_2=="P"  & valueseh=="" 
* two values blank  - C00011656 & C00011655
assert valueseh!=""  if !strpos(sample, "C0001165")
drop valueseh_*

*seu
replace valueseu="A" if valueseu_1=="A" & valueseu_2=="A" &  valueseu=="" 
replace valueseu="P" if valueseu_1=="P" | valueseu_2=="P"  & valueseu=="" 
* two values blank  - C00011656 & C00011655
assert valueseu!=""  if !strpos(sample, "C0001165")
drop valueseu_*

reshape long
rename value cleanvalue

merge 1:1 set sample method site using temp, update
noi assert value == cleanvalue if _merge==3
noi tab site if _merge ==1
replace value = cleanvalue if value==""
noi tab site if strpos(site, "_")
drop if strpos(site, "_")
noi tab site if _merge==2 
drop _merge cleanvalue
 

* counting

keep value sample site method set
reshape wide value , i(sample site) j(method) string
noi display _N " sample sites"
bysort sample: gen count1=1 if _n==1
summ count1
noi di r(sum) " samples"
bysort site: gen count2=1 if _n==1
summ count2
noi di r(sum) " sites"
drop count*

save temp2, replace


**************************************************
* SITE TYPES
***************************************************
use temp2, clear

gen type = "Phenotype" if method=="gold"

replace type = "mlst"  if inlist(site, "mlst")

replace type = "virulence" if inlist(site, "arca","arcb", "arcc", "arcd", "sasx")
replace type = "virulence" if inlist(site, "sea", "sea1", "sea2", "seb", "sec", "sed" )
replace type = "virulence" if inlist(site, "seh2", "sei", "sej", "selr", "sep","seu", "seu1", "seu2")
replace type = "virulence" if inlist(site, "eta", "etb", "etd", "chp", "sak", "scn", "lukpvf" )
replace type = "virulence" if inlist(site, "lukm", "lukmf", "lukpvs", "tsst1")
replace type = "virulence" if inlist(site, "see", "seg", "seh", "seh1" )


replace type = "ccr" if strpos(site, "ccr")==1

replace type = "Chromosomal Resistance" if inlist(site, "gyra", "grla", "dfrb", "fusa", "rpob")


replace type = "Aquired Resistance"  if inlist(site, "blaz", "meca", "mecc", "aac6aph2" )
replace type = "Aquired Resistance"  if inlist(site,"aaddaph4ia","aadeant6ia", "ant9ia", "ant9ib","aph2ic", "apha3aph3iii" )
replace type = "Aquired Resistance"  if inlist(site, "str", "fusb", "fusc", "vana", "vanb", "vanc")
replace type = "Aquired Resistance"  if inlist(site,"lnua", "lnub", "isab", "vgaa", "vgab", "vgba","erma1", "erma2", "erma")
replace type = "Aquired Resistance"  if inlist(site, "ermb", "ermc", "ermy", "ermt", "msra", "mphc", "mupa", "mupb")
replace type = "Aquired Resistance"  if inlist(site, "sat4", "tetk", "tetl", "tetm", "teto", "dfra", "dfrc", "dfrd")
replace type = "Aquired Resistance"  if inlist(site, "dfrg", "dfrk", "cfr", "cat", "qaca", "qacb", "qaccsmr")


replace type ="Phenotype" if inlist(site, "ciprofloxacin", "fusidicacid", "clindamycin", "erythromycin")
replace type ="Phenotype" if inlist(site, "gentamicin", "mupirocin", "penicillin","rifampicin")
replace type ="Phenotype" if inlist(site,"tetracycline", "trimethoprim", "vancomycin", "methicillin")

* Summary 
noi display _N " sample sites"
noi tab type
bysort sample: gen count1=1 if _n==1
summ count1
noi di r(sum) " samples"
bysort site: gen count2=1 if _n==1
summ count2
noi di r(sum) " sites"
noi tab type if count2==1
drop count*



**********************
*2) LOOK AT MISSING DATA BETWEEN METHODS
****************************************************
noi di _n(5) _dup(80) "=" _n " 2 Missing data between methods" _n _dup(80) "="

* how many values are missing
use temp2, clear
drop comments
reshape wide value, i(sample site) j(method) string

* create absence markers
gen  zbyte=( valuezam!="")
gen  tbyte=( valuetypewriter !="")
gen  gbyte=( valuegenefinder   !="")
gen new= string( tbyte ) + string(zbyte) + string(gbyte)
replace new ="Phenotype" if type=="Phenotype"
noi di "sample-sites missing values - key (typewriter mykrobe genefinder)"
table new


*000
summ zbyte if new=="000"
noi di r(N) " not called by any method; people missing from T/G and called blank by mycrobe - see 010 bonus people below"


*010
summ zbyte if new== "010"
noi di r(N) " Called only by mykrobe"
* vgaalc
noi tab site new if inlist(site, "vgaalc"), m
summ zbyte if inlist(site, "vgaalc")
noi di r(N) " vgaalc only supplied by mykrobe, drop from database"
noi drop if inlist(site, "vgaalc")
noi di _N " sample sites remaining" 
* luk
noi tab site new if strpos(site, "luk"), m
summ zbyte if inlist(site, "luk")
noi di r(N) " luk only supplied by mykrobe, drop from database"
noi drop if inlist(site, "luk")
noi di _N " sample sites remaining" 

* samples not in other sets
vallist sample if new=="010", local(list)
noi tab sample new if strpos("`list'", sample)
summ zbyte if strpos("`list'", sample)
noi di r(N) " samples are only processed by mykrobe - drop from comparision set"
drop if strpos("`list'", sample)

*001
summ zbyte if new =="001"
noi di r(N) "Called only by typewriter"
noi di "sites in sample-sites only supplied by typewriter"
noi tab site if new =="001"
* allelicprofile - drop
noi tab site new if inlist(site, "allelicprofile"), m
summ zbyte if inlist(site, "allelicprofile")
noi di r(N) " allelicprofile only supplied by typewriter, drop from database"
noi drop if inlist(site, "allelicprofile")
noi di _N " sample sites remaining" 
assert new!="001"



* 101
summ zbyte if new =="101"
noi di r(N) " no value in mykrobe"
noi tab site if new=="101"
noi tab valuezam if site=="fusa"
noi tab valuegenefinder valuetypewriter if site=="fusa" & valuezam==""
noi di " problem is that mykrobe is not reporting absence of mutation as other than blank"
vallist site if new=="101" & site!="mlst", local(list)
summ zbyte if strpos("`list'", site) & valuezam==""
noi di "solution: replace " r(N) " blanks in these sites in mykroke with wt -  not mlst which is not a mutation site"
noi replace valuezam="wt" if  strpos("`list'", site) & valuezam==""

* regen new
drop new
drop *byte
gen  zbyte=( valuezam!="")
gen  tbyte=( valuetypewriter !="")
gen  gbyte=( valuegenefinder   !="")
gen new= string( tbyte ) + string(zbyte) + string(gbyte)
noi tab new
* mslt
noi assert site=="mlst" if new=="101" 
assert new =="111" if type!="Phenotype" & site!="mlst"
noi di "all remaining sample sites are either complete accross all three methods, phenotype results, or mlst (mykrobe absent)"

* count how many sample sites left

noi display _N " sample sites"
bysort sample: gen count1=1 if _n==1
summ count1
noi di r(sum) " samples"
bysort site: gen count2=1 if _n==1
summ count2
noi di r(sum) " sites"
bysort site: gen count3=1 if _n==1 & site!="mlst" & type!="Phenotype"
summ count3
noi di r(sum) " sites not including Phenotype results and mlst"
drop count* new *byte


 
save pipeline_noncontam2, replace
*************************************



********************************
* use D-test to fix clindamycin 
replace cleanvalue=upper(value) if cleanvalue==""
reshape wide value clean comments type, i(set sample method) j(site) string
replace  cleanvalueclindamycin= "R" if  cleanvalueclindamycin=="S" &  cleanvalueclinddtest == "POS" 
drop *clinddtest
reshape long
*drop if cleanvalue==""


save pipeline_noncontam3, replace




* then add which go with which antibiotic (based on new sheet by clare)

gen antibiotic =""
replace antibiotic = "ciprofloxacin" if inlist(site, "gyra","gryb", "grla", "grlb")
replace antibiotic = "erythromycinclindamycin" if inlist(site, "erma", "ermb", "ermc", "ermt", "ermy")
replace antibiotic = "erythromycin" if inlist(site, "msra", "mphc")
replace antibiotic = "fusidic acid" if inlist(site, "fusb", "fusc", "fusa", "far")
replace antibiotic = "gentamicin" if inlist(site,"aac6aph2", "aph2ic")
replace antibiotic = "methicillinpenicillin" if inlist(site,"meca", "mecc")
replace antibiotic = "mupirocin" if inlist(site, "mupa", "mupb")
replace antibiotic = "penicillin" if site=="blaz"
replace antibiotic = "rifampicin" if site=="rpob"
replace antibiotic = "tetracycline" if inlist(site, "tetk", "tetl", "tetm", "teto")
replace antibiotic = "trimethoprim" if inlist(site, "dfra", "dfrc", "dfrb", "dfrb","dfrg", "dfrk")
replace antibiotic = "vancomycin" if inlist(site, "vana", "vanb", "vanc")



sort set sample method site
save pipeline_noncontam4, replace



*clean up non binary values
use pipeline_noncontam4, clear
gen orgvalue = value
replace cleanvalue = upper(value) if cleanvalue==""
replace cleanvalue="" if inlist(value, "nd", "ND", "NT", "Low cov", "mixed")
replace cleanvalue="" if inlist(cleanvalue, "NOT DONE", "N/A", "NF")

replace cleanvalue="" if type=="profile"
replace cleanvalue="" if type=="mlst"

*I/B are equivalent = somewhere between sensitive and resistant
replace cleanvalue="I" if value=="B"


replace cleanvalue="P" if inlist(value, "Pi", "Pu")
replace cleanvalue="A" if value=="Pc"

replace value = subinstr(value, "=>", "-",.)
replace cleanvalue = subinstr(cleanvalue, "=>", "-",.)
replace cleanvalue = subinstr(cleanvalue, " ", "",.)

*swap to binary the point mutations
replace cleanvalue = "A" if cleanvalue=="WT"
replace cleanvalue = "A" if cleanvalue=="" & method=="zam" & strpos(type, "Chromo")

* point mutations based on clare's paper
*gyra
replace cleanvalue = "P" if cleanvalue =="84:S-L" & strpos(site, "gyra")
*grla 
replace cleanvalue = "P" if strpos(cleanvalue, "80:S-F") & strpos(site, "grla")
replace cleanvalue = "P" if strpos(cleanvalue, "80:S-Y") & strpos(site, "grla")
replace cleanvalue = "P" if strpos(cleanvalue, "84:E-G") & strpos(site, "grla")

* fusa
replace cleanvalue = "P" if strpos(cleanvalue, "114:P-H") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "461:L-K") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "461:L-S") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "404:P-L") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "90:V-I") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "457:H-Y") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "67:A-T") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "406:P-L") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "464:R-S") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "464:R-H") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "457:H-Q") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "453:M-I") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "452:G-C") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "452:G-S") & strpos(site, "fusa")
replace cleanvalue = "P" if strpos(cleanvalue, "404:P-Q") & strpos(site, "fusa")

*rpob
replace cleanvalue = "P" if strpos(cleanvalue, "464:S-P") & strpos(site, "rpob")
replace cleanvalue = "P" if strpos(cleanvalue, "468:Q-R") & strpos(site, "rpob")
replace cleanvalue = "P" if strpos(cleanvalue, "471:D-Y") & strpos(site, "rpob")
replace cleanvalue = "P" if strpos(cleanvalue, "477:A-D") & strpos(site, "rpob")
replace cleanvalue = "P" if strpos(cleanvalue, "481:H-N") & strpos(site, "rpob")
replace cleanvalue = "P" if strpos(cleanvalue, "481:H-Y") & strpos(site, "rpob")

*dfrb

replace cleanvalue = "P" if strpos(cleanvalue, "150:H-R") & strpos(site, "dfrb")
replace cleanvalue = "P" if strpos(cleanvalue, "21:L-V") & strpos(site, "dfrb")
replace cleanvalue = "P" if strpos(cleanvalue, "31:H-N") & strpos(site, "dfrb")
replace cleanvalue = "P" if strpos(cleanvalue, "41:L-F") & strpos(site, "dfrb")
replace cleanvalue = "P" if strpos(cleanvalue, "99:F-S") & strpos(site, "dfrb")
replace cleanvalue = "P" if strpos(cleanvalue, "99:F-Y") & strpos(site, "dfrb")



* replace all other point mutations as A

replace cleanvalue = "A" if cleanvalue!="P"& strpos(type, "Chromo")


* swap blanks to dashes to make easier to see
replace cleanvalue="-" if cleanvalue==""


save pipeline_clean, replace


