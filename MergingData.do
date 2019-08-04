***STATA code to merge surveys on phthalates, phenols and parabens, and pregnancy data 2005-2014
***Written with Stata 13.1


*** SURVEY D, 2005-2006
*import pregnancy data
clear all
import sasxport "UCPREG_D.XPT"
tab urxpreg
keep if urxpreg==1
qui summ seqn
assert r(N)==367 
save UCPREG_D, replace
*import phthalates data
clear all
import sasxport "PHTHTE_D.XPT"
summ
qui summ seqn
assert r(N)==2638 
save PHTHTE_D, replace
*import phenols and parabens data, merge with pthalates and pregnancy data
clear all
import sasxport "EPH_D.XPT"
merge 1:1 seqn wtsb2yr urxucr using "PHTHTE_D.dta"
drop if _merge==1 /* two obsns have diff urxucr, keep either one as both have missing exp data*/
summ
qui summ seqn
assert r(N)==2638 
drop _merge
merge 1:1 seqn using "UCPREG_D.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EHP_UCPREG_D, replace


*** SURVEY E, 2007-2008
*import pregnancy data
clear all
import sasxport "UCPREG_E.XPT"
tab urxpreg
keep if urxpreg==1
qui summ seqn
assert r(N)==57 
save UCPREG_E, replace
* import phthalates data
clear all
import sasxport "PHTHTE_E.XPT"
summ
qui summ seqn
assert r(N)==2718 
save PHTHTE_E, replace
*import phenols and parabens data, merge with pthalates and pregnancy data
clear all
import sasxport "EPH_E.XPT"
qui summ seqn
assert r(N)==2718 
merge 1:1 seqn wtsb2yr urxucr using "PHTHTE_E.dta"
summ
drop _merge
merge 1:1 seqn using "UCPREG_E.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EHP_UCPREG_E, replace


*** SURVEY F, 2009-2010
*import pregnancy data
clear all
import sasxport "UCPREG_F.XPT"
tab urxpreg
keep if urxpreg==1
qui summ seqn
assert r(N)==68 
save UCPREG_F, replace
* import phthalates data
clear all
import sasxport "PHTHTE_F.XPT"
summ
qui summ seqn
assert r(N)==2819 
save PHTHTE_F, replace
*import phenols and parabens data, merge with pthalates and pregnancy data
clear all
import sasxport "EPH_F.XPT"
qui summ seqn
assert r(N)==2819
merge 1:1 seqn wtsb2yr urxucr using "PHTHTE_F.dta"
summ
drop _merge
merge 1:1 seqn using "UCPREG_F.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EHP_UCPREG_F, replace


*** SURVEY G, 2011-2012
*import pregnancy data 
clear all
import sasxport "UCPREG_G.XPT"
tab urxpreg
keep if urxpreg==1
qui summ seqn
assert r(N)==55 
save UCPREG_G, replace
* import phthalates data
clear all
import sasxport "PHTHTE_G.XPT"
summ
qui summ seqn
assert r(N)==2594 
save PHTHTE_G, replace
*import phenols and parabens data, merge with pthalates and pregnancy data
clear all
import sasxport "EPH_G.XPT"
qui summ seqn
assert r(N)==2594
merge 1:1 seqn wtsb2yr urxucr using "PHTHTE_G.dta"
summ
drop _merge
merge 1:1 seqn using "UCPREG_G.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EHP_UCPREG_G, replace


*** SURVEY H, 2013-2014
*import pregnancy data
clear all
import sasxport "UCPREG_H.XPT"
tab urxpreg
keep if urxpreg==1
qui summ seqn
assert r(N)==60 
save UCPREG_H, replace
* import phthalates data
clear all
import sasxport "PHTHTE_H.XPT"
summ
qui summ seqn
assert r(N)==2777
save PHTHTE_H, replace
*import phenols and parabens data, merge with pthalates and pregnancy data
clear all
import sasxport "EPHPP_H.XPT"
qui summ seqn
assert r(N)==2777
merge 1:1 seqn urxucr using "PHTHTE_H.dta"
summ
drop _merge
merge 1:1 seqn using "UCPREG_H.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EPHPP_UCPREG_H, replace


***LOAD MERGED DATA FROM SURVEYS D-H 2005-2014
clear all
use NHANESdata_PHTHTE_EHP_UCPREG_D
gen survey="D"
summ 
append using NHANESdata_PHTHTE_EHP_UCPREG_E
replace survey="E" if survey==""
summ
append using NHANESdata_PHTHTE_EHP_UCPREG_F
replace survey="F" if survey==""
summ
append using NHANESdata_PHTHTE_EHP_UCPREG_G
replace survey="G" if survey==""
summ 
append using NHANESdata_PHTHTE_EPHPP_UCPREG_H
replace survey="H" if survey==""
summ

*create indicator for below max LOD across all surveys
gen urdcoplcmax=urdcoplc
replace urdcoplcmax=1 if urxcop<0.7 & urdcoplc==0
gen urdcnplcmax=urdcnplc
replace urdcnplcmax=1 if urxcnp<0.6 & urdcnplc==0
gen urdecplcmax=urdecplc
replace urdecplcmax=1 if urxecp<0.6 & urdecplc==0
gen urdmbplcmax=urdmbplc
replace urdmbplcmax=1 if urxmbp<0.6 & urdmbplc==0
gen urdmc1lcmax=urdmc1lc
replace urdmc1lcmax=1 if urxmc1<0.4 & urdmc1lc==0
gen urdmeplcmax=urdmeplc
replace urdmeplcmax=1 if urxmep<1.2 & urdmeplc==0
gen urdmhhlcmax=urdmhhlc
replace urdmhhlcmax=1 if urxmhh<0.7 & urdmhhlc==0
gen urdmhplcmax=urdmhplc 
replace urdmhplcmax=1 if urxmhp<1.2 & urdmhplc==0
gen urdmiblcmax=urdmiblc
replace urdmiblcmax=1 if urxmib<0.8 & urdmiblc==0
gen urdmnplcmax=urdmnplc
replace urdmnplcmax=1 if urxmnp<1.232 & urdmnplc==0
gen urdmohlcmax=urdmohlc
replace urdmohlcmax=1 if urxmoh<0.7 & urdmohlc==0
gen urdmzplcmax=urdmzplc
replace urdmzplcmax=1 if urxmzp<0.3 & urdmzplc==0
gen urdbphlcmax=urdbphlc
replace urdbphlcmax=1 if urxbph<0.4 & urdbphlc==0
gen urdtrslcmax=urdtrslc 
replace urdtrslcmax=1 if urxtrs<2.3 & urdtrslc==0
gen urdbp3lcmax=urdbp3lc
replace urdbp3lcmax=1 if urxbp3<0.4 & urdbp3lc==0
gen urdppblcmax=urdppblc
replace urdppblcmax=1 if urxppb<0.2 & urdppblc==0
gen urdbuplcmax=urdbuplc
replace urdbuplcmax=1 if urxbup<0.2 & urdbuplc==0
gen urdepblcmax=urdepblc 
replace urdepblcmax=1 if urxepb<1 & urdepblc==0
gen urdmpblcmax=urdmpblc
replace urdmpblcmax=1 if urxmpb<1 & urdmpblc==0

*replace below maxLOD values with maxLOD/sqrt(2) 
gen urxcop_mlod=urxcop
replace urxcop_mlod=0.7/sqrt(2) if urdcoplcmax==1
gen urxcnp_mlod=urxcnp
replace urxcnp_mlod=0.6/sqrt(2) if urdcnplcmax==1 
gen urxecp_mlod=urxecp
replace urxecp_mlod=0.6/sqrt(2) if urdecplcmax==1 
gen urxmbp_mlod=urxmbp
replace urxmbp_mlod=0.6/sqrt(2) if urdmbplcmax==1 
gen urxmc1_mlod=urxmc1
replace urxmc1_mlod=0.4/sqrt(2) if urdmc1lcmax==1 
gen urxmep_mlod=urxmep
replace urxmep_mlod=1.2/sqrt(2) if urdmeplcmax==1 
gen urxmhh_mlod=urxmhh
replace urxmhh_mlod=0.7/sqrt(2) if urdmhhlcmax==1 
gen urxmhp_mlod=urxmhp 
replace urxmhp_mlod=1.2/sqrt(2) if urdmhplcmax==1 
gen urxmib_mlod=urxmib
replace urxmib_mlod=0.8/sqrt(2) if urdmiblcmax==1 
gen urxmnp_mlod=urxmnp
replace urxmnp_mlod=1.232/sqrt(2) if urdmnplcmax==1 
gen urxmoh_mlod=urxmoh
replace urxmoh_mlod=0.7/sqrt(2) if urdmohlcmax==1 
gen urxmzp_mlod=urxmzp
replace urxmzp_mlod=0.3/sqrt(2) if urdmzplcmax==1 
gen urxbph_mlod=urxbph
replace urxbph_mlod=0.4/sqrt(2) if urdbphlcmax==1 
gen urxtrs_mlod=urxtrs 
replace urxtrs_mlod=2.3/sqrt(2) if urdtrslcmax==1 
gen urxbp3_mlod=urxbp3
replace urxbp3_mlod=0.4/sqrt(2) if urdbp3lcmax==1 
gen urxppb_mlod=urxppb
replace urxppb_mlod=0.2/sqrt(2) if urdppblcmax==1 
gen urxbup_mlod=urxbup
replace urxbup_mlod=0.2/sqrt(2) if urdbuplcmax==1 
gen urxepb_mlod=urxepb 
replace urxepb_mlod=1/sqrt(2) if urdepblcmax==1 
gen urxmpb_mlod=urxmpb
replace urxmpb_mlod=1/sqrt(2) if urdmpblcmax==1 

*rename BPH to BPA, common name
rename urdbphlcmax urdbpalcmax
rename urxbph_mlod urxbpa_mlod

*create creatinine-adjusted exposure vars
gen cop_ucr=urxcop/urxucr
gen cnp_ucr=urxcnp/urxucr
gen ecp_ucr=urxecp/urxucr 
gen mbp_ucr=urxmbp/urxucr 
gen mc1_ucr=urxmc1/urxucr 
gen mep_ucr=urxmep/urxucr 
gen mhh_ucr=urxmhh/urxucr 
gen mhp_ucr=urxmhp/urxucr 
gen mib_ucr=urxmib/urxucr 
gen mnp_ucr=urxmnp/urxucr 
gen moh_ucr=urxmoh/urxucr 
gen mzp_ucr=urxmzp/urxucr 
gen bpa_ucr=urxbph/urxucr 
gen trs_ucr=urxtrs/urxucr 
gen bp3_ucr=urxbp3/urxucr 
gen ppb_ucr=urxppb/urxucr 
gen bup_ucr=urxbup/urxucr 
gen epb_ucr=urxepb/urxucr 
gen mpb_ucr=urxmpb/urxucr 

*create creatinine-adjusted exposure vars, using max LOD
gen cop_mlod_ucr=urxcop_mlod/urxucr
gen cnp_mlod_ucr=urxcnp_mlod/urxucr
gen ecp_mlod_ucr=urxecp_mlod/urxucr 
gen mbp_mlod_ucr=urxmbp_mlod/urxucr 
gen mc1_mlod_ucr=urxmc1_mlod/urxucr 
gen mep_mlod_ucr=urxmep_mlod/urxucr 
gen mhh_mlod_ucr=urxmhh_mlod/urxucr 
gen mhp_mlod_ucr=urxmhp_mlod/urxucr 
gen mib_mlod_ucr=urxmib_mlod/urxucr 
gen mnp_mlod_ucr=urxmnp_mlod/urxucr 
gen moh_mlod_ucr=urxmoh_mlod/urxucr 
gen mzp_mlod_ucr=urxmzp_mlod/urxucr 
gen bpa_mlod_ucr=urxbpa_mlod/urxucr 
gen trs_mlod_ucr=urxtrs_mlod/urxucr 
gen bp3_mlod_ucr=urxbp3_mlod/urxucr 
gen ppb_mlod_ucr=urxppb_mlod/urxucr 
gen bup_mlod_ucr=urxbup_mlod/urxucr 
gen epb_mlod_ucr=urxepb_mlod/urxucr 
gen mpb_mlod_ucr=urxmpb_mlod/urxucr 

*drop chems not measured in all 5 surveys
drop *mcp* *mnm* *mop* *4to*

*drop 4 observations with all exposure measurements missing
drop if urxcop==. & urxcnp==. & urxecp==. & urxmbp==. & urxmc1==. & urxmep==. & urxmhh==. & urxmhp==. & urxmib==. & urxmnp==. & urxmoh==. & urxmzp==. & urxbph==. & urxtrs==. & urxbp3==. & urxppb==. & urxbup==. & urxepb==. & urxmpb==. 

*label new vars
label var cop_ucr "Urinary Mono(carboxyoctyl) phthalate (ug/g creatinine)"
label var cnp_ucr "Urinary Mono(carboxynonyl) phthalate (ug/g creatinine)"
label var ecp_ucr "Urinary Mono-2-ethyl-5-carboxypentyl phthalate (ug/g creatinine)"
label var mbp_ucr "Urinary Mono-n-butyl phthalate (ug/g creatinine)"
label var mc1_ucr "Urinary Mono-(3-carboxypropyl) phthalate (ug/g creatinine)"
label var mep_ucr "Urinary Mono-ethyl phthalate (ug/g creatinine)"
label var mhh_ucr "Urinary Mono-(2-ethyl-5-hydroxyhexyl) phthalate (ug/g creatinine)"
label var mhp_ucr "Urinary Mono-2-ethyl-5-carboxypentyl phthalate (ug/g creatinine)"
label var mib_ucr "Urinary Mono-isobutyl phthalate (ug/g creatinine)"
label var mnp_ucr "Urinary Mono-isononyl phthalate (ug/g creatinine)"
label var moh_ucr "Urinary Mono-(2-ethyl-5-oxohexyl) phthalate (ug/g creatinine)"
label var mzp_ucr "Urinary Mono-benzyl phthalate (ug/g creatinine)"
label var bpa_ucr "Urinary Bisphenol A (ug/g creatinine)"
label var trs_ucr "Urinary Triclosan (ug/g creatinine)" 
label var bp3_ucr "Urinary Benzophenone-3 (ug/g creatinine)"
label var ppb_ucr "Urinary Propyl paraben (ug/g creatinine)"
label var bup_ucr "Urinary Butyl paraben (ug/g creatinine)"
label var epb_ucr "Urinary Ethyl paraben (ug/g creatinine)" 
label var mpb_ucr "Urinary Methyl paraben (ug/g creatinine)" 

*label new vars
label var cop_mlod_ucr "Urinary Mono(carboxyoctyl) phthalate (ug/g creatinine)"
label var cnp_mlod_ucr "Urinary Mono(carboxynonyl) phthalate (ug/g creatinine)"
label var ecp_mlod_ucr "Urinary Mono-2-ethyl-5-carboxypentyl phthalate (ug/g creatinine)"
label var mbp_mlod_ucr "Urinary Mono-n-butyl phthalate (ug/g creatinine)"
label var mc1_mlod_ucr "Urinary Mono-(3-carboxypropyl) phthalate (ug/g creatinine)"
label var mep_mlod_ucr "Urinary Mono-ethyl phthalate (ug/g creatinine)"
label var mhh_mlod_ucr "Urinary Mono-(2-ethyl-5-hydroxyhexyl) phthalate (ug/g creatinine)"
label var mhp_mlod_ucr "Urinary Mono-2-ethyl-5-carboxypentyl phthalate (ug/g creatinine)"
label var mib_mlod_ucr "Urinary Mono-isobutyl phthalate (ug/g creatinine)"
label var mnp_mlod_ucr "Urinary Mono-isononyl phthalate (ug/g creatinine)"
label var moh_mlod_ucr "Urinary Mono-(2-ethyl-5-oxohexyl) phthalate (ug/g creatinine)"
label var mzp_mlod_ucr "Urinary Mono-benzyl phthalate (ug/g creatinine)"
label var bpa_mlod_ucr "Urinary Bisphenol A (ug/g creatinine)"
label var trs_mlod_ucr "Urinary Triclosan (ug/g creatinine)" 
label var bp3_mlod_ucr "Urinary Benzophenone-3 (ug/g creatinine)"
label var ppb_mlod_ucr "Urinary Propyl paraben (ug/g creatinine)"
label var bup_mlod_ucr "Urinary Butyl paraben (ug/g creatinine)"
label var epb_mlod_ucr "Urinary Ethyl paraben (ug/g creatinine)" 
label var mpb_mlod_ucr "Urinary Methyl paraben (ug/g creatinine)"

*save
save NHANES_ExposureData, replace
export delimited NHANES_ExposureData, replace

*keep only exposures with greater than 90% detection rate
summ urdcoplcmax urdcnplcmax urdecplcmax urdmbplcmax urdmc1lcmax urdmeplcmax urdmhhlcmax urdmhplcmax urdmiblcmax urdmnplcmax urdmohlcmax urdmzplcmax urdbpalcmax urdtrslcmax urdbp3lcmax urdppblcmax urdbuplcmax urdepblcmax urdmpblcmax 
ds *max
local maxvarlist `"`r(varlist)'"'
foreach var in `maxvarlist' {
	summ `var', meanonly
	if r(mean)<0.1 {
		local varabb=substr("`var'",4,3)
		local keeplist "`keeplist' `varabb'_mlod_ucr"
	}
}
keep `keeplist'

*rename to only 3 letter abbreviation
qui ds
local varlist `"`r(varlist)'"'
foreach var in `varlist' {	
	local varabb=upper(substr("`var'",1,3))
	rename `var' `varabb'
}

*save
save NHANES_ExposureData_12, replace
export delimited NHANES_ExposureData_12, replace


