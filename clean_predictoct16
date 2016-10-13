* clean pipeline data

set li 130

cap log close
log using clean.log, replace
noi di "Run by AMM on $S_DATE $S_TIME"
cd E:\users\amy.mason\Pipeline_27_07_2016\Datasets

noi di _n(5) _dup(80) "=" _n " 1 Combine all predictive sets" _n _dup(80) "="


****** combine the 3 sets
noi di "combine"
use pipeline_data_gf, clear
append using pipeline_data_tw.dta, force
append using pipeline_data_z2.dta, force
gsort sample - set
by sample: replace set=set[_n-1] if _n>1 & set==""
by sample: assert set==set[_n-1] if _n>1
replace site = lower(site)
noi di _N " results"
bysort site: gen tab=1 if _n==1
summ tab 
noi di "on " r(sum) " sites"
tab site
drop tab
bysort sample: gen tab=1 if _n==1
summ tab
noi di "from " r(sum) " samples"
tab set if tab==1
drop tab

save pipeline_predict_raw, replace

***********************************************************
* 1) Remove samples that are contaminated for all methods
*********************************************************
use pipeline_predict_raw, clear

noi di _n(5) _dup(80) "=" _n " 1 Remove samples that were contaminated for all methods" _n _dup(80) "="

use pipeline_predict_raw, clear
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
assert _merge!=4
assert _merge!=5
noi assert value == cleanvalue if _merge==3
noi tab site if _merge ==1
replace value = cleanvalue if value==""
noi tab site if strpos(site, "_")
drop if strpos(site, "_")
noi tab site if _merge==2 
drop _merge cleanvalue comment
 

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

gen type = "mlst"  if inlist(site, "mlst")

replace type = "virulence" if inlist(site, "arca","arcb", "arcc", "arcd", "sasx")
replace type = "virulence" if inlist(site, "sea", "sea1", "sea2", "seb", "sec", "sed" )
replace type = "virulence" if inlist(site, "seh2", "sei", "sej", "selr", "sep","seu", "seu1", "seu2")
replace type = "virulence" if inlist(site, "eta", "etb", "etd", "chp", "sak", "scn", "lukpvf" )
replace type = "virulence" if inlist(site, "luk", "lukm", "lukmf", "lukpvs", "tsst1")
replace type = "virulence" if inlist(site, "see", "seg", "seh", "seh1" )


replace type = "ccr" if strpos(site, "ccr")==1

replace type = "Chromosomal Resistance" if inlist(site, "iles", "gyra", "grla", "dfrb", "fusa", "rpob")


replace type = "Aquired Resistance"  if inlist(site, "blaz", "meca", "mecc", "aac6aph2" )
replace type = "Aquired Resistance"  if inlist(site,"aaddaph4ia","aadeant6ia", "ant9ia", "ant9ib","aph2ic", "apha3aph3iii" )
replace type = "Aquired Resistance"  if inlist(site, "str", "fusb", "fusc", "vana", "vanb", "vanc")
replace type = "Aquired Resistance"  if inlist(site,"lnua", "lnub", "isab", "vgaa", "vgab", "vgba","erma1", "erma2", "erma")
replace type = "Aquired Resistance"  if inlist(site, "ermb", "ermc", "ermy", "ermt", "msra", "mphc", "mupa", "mupb")
replace type = "Aquired Resistance"  if inlist(site, "sat4", "tetk", "tetl", "tetm", "teto", "dfra", "dfrc", "dfrd")
replace type = "Aquired Resistance"  if inlist(site, "dfrg", "dfrk", "cfr", "cat", "qaca", "qacb", "qaccsmr")



* Summary 
noi display _N " sample sites"
noi tab type
bysort sample: gen count1=1 if _n==1
summ count1
noi di r(sum) " samples"
bysort site: gen count2=1 if _n==1
summ count2
noi di r(sum) " sites"
noi tab type if count2==1, m
noi tab site if type==""
drop count*

save temp3, replace


**********************
*2) LOOK AT MISSING DATA BETWEEN METHODS
****************************************************
noi di _n(5) _dup(80) "=" _n " 2 Missing data between methods" _n _dup(80) "="

* how many values are missing
use temp3, clear

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
* ar
noi tab site new if strpos(site, "ar"), m
summ zbyte if inlist(site, "ar")
noi di r(N) " ar only supplied by mykrobe, drop from database"
noi drop if inlist(site, "ar")
noi di _N " sample sites remaining" 


* samples not in other sets
vallist sample if new=="010", local(list)
noi tab sample new if strpos("`list'", sample)
bysort sample: gen count =1 if _n==1
summ zbyte if strpos("`list'", sample) & count==1
noi di r(N) " samples are only processed by mykrobe - drop from comparision set"
summ zbyte if strpos("`list'", sample) 
noi di r(N) " sample sites affected"
drop if strpos("`list'", sample)
drop count 

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
noi di r(N) " no value in mykrobe, eg fusa"
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
assert new =="111" if (type!="Phenotype" & site!="mlst")
noi di "all remaining sample sites are either complete accross all three methods or mlst (mykrobe absent)"

* save mlst results and drop
preserve
keep if site =="mlst"
gen ident =  strpos(valuetypewriter, valuegene) >0
replace ident = 1 if valueg=="Novel" & valuet=="NF"
summ ident 
noi di "mlst agreement between typewriter and genefinder in " r(sum) " cases out of " r(N)
noi di "value missing for all mykrobe"
noi di "rest"
noi tab valueg valuet if ident==0
keep sample valueg valuet ident
save mlst, replace
restore
gen marker =1 if site=="mlst"
summ marker
noi di "dropping " r(sum) " mlst sites"
noi drop if marker==1
drop marker


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


 
save temp4, replace
*************************************

* then add which go with which antibiotic (based on new sheet by clare)
**********************************************************

noi di _n(5) _dup(80) "=" _n " 4 Match to antibiotic prediction" _n _dup(80) "="


*clean up non binary values
use temp4, clear
reshape long


 noi di "missing values"
replace value = upper(value)
gen ambi = 1 if inlist(value, "ND", "NT", "LOW COV", "MIXED")
replace ambi=1 if inlist(value, "NOT DONE", "N/A", "NF", "NA", "-", "X")
summ ambi if ambi ==1
noi di r(N) " values have N/A or equivalent values, set these all to A"
noi list site value if ambi==1
noi replace value = "A" if ambi==1
tab site if ambi ==1
drop ambi

noi di "all values equal A or P, unless chromosonal"
replace value=upper(value)

assert inlist(value, "A", "P") if !strpos(type, "Chro")

***************
gen count = 1 if strpos(value, "=>")
summ count
noi di r(N) " replace => with - to create uniform recording of mutations"
noi replace value = subinstr(value, "=>", "-",.)
drop count 


gen count = 1 if strpos(value," ")
summ count
noi di r(N) " remove all spaces in value reports"
noi replace value = subinstr(value, " ", "",.)
drop count 

* fix new mykrobe output
gen start=strpos(value,"_")
gen end= strpos(value, "-")
replace value= substr(value,start+2,end-start-3)+":"+substr(value,start+1,1)+"-" +substr(value,end-1,1) if site==lower(substr(value,1,4)) & strpos(type, "Chr") & method=="zam"
drop start end

*swap to binary the point mutations
* clarify absenses
gen count = 1  if value=="WT" & strpos(type, "Chromo")
summ count if count==1 
noi di r(N) "replace WT with A in typewriter/ genewriter/ zam  reports on chromosonal mutations"
noi replace value = "A" if count==1 
drop count

* clarify point mutations based on clare's paper
noi di "mutations: if count!=1 then no considered positive as not in panal on clare's paper"
*gyra
noi di "gyra mutations"
gen count =1 if value =="84:S-L" & strpos(site, "gyra")
noi tab value count if site=="gyra" & strpos(value, "-") , m
replace value = "P" if count ==1
noi tab value count if site=="gyra" , m
drop count 
*grla
noi di "grla mutations"
gen count =1 if strpos(value, "80:S-F") & strpos(site, "grla")
replace count=1 if strpos(value, "80:S-Y") & strpos(site, "grla")
replace count=1  if strpos(value, "84:E-G") & strpos(site, "grla")
noi tab value count if site=="grla" & strpos(value, "-") , m
replace value = "P" if count ==1
noi tab value count if site=="grla" , m
drop count 


* fusa
noi di "fusa mutations"

gen count =1 if strpos(value, "114:P-H") & strpos(site, "fusa")
replace count =1  if strpos(value, "461:L-K") & strpos(site, "fusa")
replace count =1  if strpos(value, "461:L-S") & strpos(site, "fusa")
replace count =1  if strpos(value, "404:P-L") & strpos(site, "fusa")
replace count =1  if strpos(value, "90:V-I") & strpos(site, "fusa")
replace count =1  if strpos(value, "457:H-Y") & strpos(site, "fusa")
replace count =1 if strpos(value, "67:A-T") & strpos(site, "fusa")
replace count =1 if strpos(value, "406:P-L") & strpos(site, "fusa")
replace count =1 if strpos(value, "464:R-S") & strpos(site, "fusa")
replace count =1  if strpos(value, "464:R-H") & strpos(site, "fusa")
replace count =1  if strpos(value, "457:H-Q") & strpos(site, "fusa")
replace count =1  if strpos(value, "453:M-I") & strpos(site, "fusa")
replace count =1  if strpos(value, "452:G-C") & strpos(site, "fusa")
replace count =1  if strpos(value, "452:G-S") & strpos(site, "fusa")
replace count =1  if strpos(value, "404:P-Q") & strpos(site, "fusa")
noi tab value count if site=="fusa" & strpos(value, "-") , m
replace value = "P" if count ==1
noi tab value count if site=="fusa" , m
drop count 


*rpob
noi di "rpob mutations"
gen count =1 if strpos(value, "464:S-P") & strpos(site, "rpob")
replace count =1   if strpos(value, "468:Q-R") & strpos(site, "rpob")
replace count =1  if strpos(value, "471:D-Y") & strpos(site, "rpob")
replace count =1  if strpos(value, "477:A-D") & strpos(site, "rpob")
replace count =1   if strpos(value, "481:H-N") & strpos(site, "rpob")
replace count =1  if strpos(value, "481:H-Y") & strpos(site, "rpob")
noi tab value count if site=="rpob" & strpos(value, "-") , m
noi replace value = "P" if count ==1
noi tab value count if site=="rpob" , m
drop count 


*dfrb
noi di "dfrb mutations"
gen count =1 if  strpos(value, "150:H-R") & strpos(site, "dfrb")
replace count =1   if  strpos(value, "21:L-V") & strpos(site, "dfrb")
replace count =1   if strpos(value, "31:H-N") & strpos(site, "dfrb")
replace count =1   if strpos(value, "41:L-F") & strpos(site, "dfrb")
replace count =1   if strpos(value, "99:F-S") & strpos(site, "dfrb")
replace count =1   if  strpos(value, "99:F-Y") & strpos(site, "dfrb")
noi tab value count if site=="dfrb" & strpos(value, "-") , m
noi replace value = "P" if count ==1
noi tab value count if site=="dfrb" , m
drop count 

* replace all other point mutations as A
noi di "replace mutations not in table with A"
gen count=1 if !inlist(value, "P", "A")& strpos(type, "Chromo")& inlist(method, "genefinder", "typewriter", "zam") & value!=""
tab value site if count==1, m
replace value="A" if count==1 & value!="NA"
drop count



* swap blanks to dashes to make easier to see
assert value!=""
replace site=lower(site)
noi bysort type: tab site value, m 

save temp5, replace

noi di _n(5) _dup(80) "=" _n " 5 Cleaned summary" _n _dup(80) "="

***************************************
* Summary 
****************************************
noi display _N " sample sites"
noi tab type
bysort sample: gen count1=1 if _n==1
summ count1
noi di r(sum) " samples"
bysort site: gen count2=1 if _n==1
summ count2
noi di r(sum) " sites"
noi tab type if count2==1, m
noi tab site if type==""
drop count*


noi di "values by type"
noi bysort type: tab value 

noi di "breakdown by site"
noi bysort type site:  tab value 

noi di "breakdown by method"
noi bysort method: tab value

save pipeline_clean_all_values_long, replace

**********
* and wide
noi di "Compare site values between methods"

reshape wide
gen ValueAll = valuet  + valuez + valueg
noi tab ValueAll , sort m

noi di "values by type"
noi bysort type: tab ValueA, sort  m 

noi di "breakdown by site"
noi bysort type site:  tab ValueA, sort m

noi tab site ValueA, m

save pipeline_clean_all_values_wide, replace

cd E:\users\amy.mason\Pipeline_27_07_2016\


