***STATA code to merge surveys on phthalates, phenols and parabens, and pregnancy data 2005-2014
***Written with Stata 13.1

***required NHANES data files:
***SURVEY D, 2005-2006: UCPREG_D.XPT, PHTHTE_D.XPT, EPH_D.XPT, PP_D.XPT
***SURVEY E, 2007-2008: UCPREG_E.XPT, PHTHTE_E.XPT, EPH_E.XPT, PP_E.XPT
***SURVEY F, 2009-2010: UCPREG_F.XPT, PHTHTE_F.XPT, EPH_F.XPT, PP_F.XPT
***SURVEY G, 2011-2012: UCPREG_G.XPT, PHTHTE_G.XPT, EPH_G.XPT, PP_G.XPT
***SURVEY H, 2013-2014: UCPREG_H.XPT, PHTHTE_H.XPT, EPHPP_H.XPT

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
*import extra phenols data
clear all
import sasxport "PP_D.XPT"
summ
qui summ seqn
assert r(N)==2638 
save PP_D, replace
*import phenols and parabens data, merge with pthalates data
clear all
import sasxport "EPH_D.XPT"
merge 1:1 seqn wtsb2yr urxucr using "PHTHTE_D.dta"
drop if _merge==1 /* two obsns have diff urxucr, keep either one as both have missing exp data*/
summ
qui summ seqn
assert r(N)==2638 
drop _merge
*merge with extra phenols data
merge 1:1 seqn wtsb2yr urxucr using "PP_D.dta"
drop if _merge==1 /* two obsns have diff urxucr, keep either one as both have missing exp data*/
summ
qui summ seqn
assert r(N)==2638 
drop _merge
*merge with pregnancy data
merge 1:1 seqn using "UCPREG_D.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EPH_PP_UCPREG_D, replace


*** SURVEY E, 2007-2008
*import pregnancy data
clear all
import sasxport "UCPREG_E.XPT"
tab urxpreg
keep if urxpreg==1
qui summ seqn
assert r(N)==57 
save UCPREG_E, replace
*import phthalates data
clear all
import sasxport "PHTHTE_E.XPT"
summ
qui summ seqn
assert r(N)==2718 
save PHTHTE_E, replace
*import extra phenols data
clear all
import sasxport "PP_E.XPT"
summ
qui summ seqn
assert r(N)==2718 
save PP_E, replace
*import phenols and parabens data, merge with pthalates data
clear all
import sasxport "EPH_E.XPT"
qui summ seqn
assert r(N)==2718 
merge 1:1 seqn wtsb2yr urxucr using "PHTHTE_E.dta"
summ
drop _merge
*merge with extra phenols data
merge 1:1 seqn wtsb2yr urxucr using "PP_E.dta"
keep if _merge==3
qui summ seqn
assert r(N)==2718 
drop _merge
*merge with pregnancy data
merge 1:1 seqn using "UCPREG_E.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EPH_PP_UCPREG_E, replace


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
*import extra phenols data
clear all
import sasxport "PP_F.XPT"
summ
qui summ seqn
assert r(N)==2819 
save PP_F, replace
*import phenols and parabens data, merge with pthalates data
clear all
import sasxport "EPH_F.XPT"
qui summ seqn
assert r(N)==2819
merge 1:1 seqn wtsb2yr urxucr using "PHTHTE_F.dta"
summ
drop _merge
*merge with extra phenols data
merge 1:1 seqn wtsb2yr urxucr using "PP_F.dta"
keep if _merge==3
qui summ seqn
assert r(N)==2819 
drop _merge
*merge with pregnancy data
merge 1:1 seqn using "UCPREG_F.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EPH_PP_UCPREG_F, replace


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
*import extra phenols data
clear all
import sasxport "PP_G.XPT.txt"
summ
qui summ seqn
assert r(N)==2594 
save PP_G, replace
*import phenols and parabens data, merge with pthalates data
clear all
import sasxport "EPH_G.XPT"
qui summ seqn
assert r(N)==2594
merge 1:1 seqn wtsa2yr urxucr using "PHTHTE_G.dta"
summ
drop _merge
*merge with extra phenols data
merge 1:1 seqn wtsa2yr urxucr using "PP_G.dta"
keep if _merge==3
qui summ seqn
assert r(N)==2594 
drop _merge
*merge with pregnancy data
merge 1:1 seqn using "UCPREG_G.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EPH_PP_UCPREG_G, replace


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
*import phenols and parabens data, merge with pthalates data
clear all
import sasxport "EPHPP_H.XPT"
qui summ seqn
assert r(N)==2777
merge 1:1 seqn urxucr using "PHTHTE_H.dta"
summ
drop _merge
*merge with pregnancy data
merge 1:1 seqn using "UCPREG_H.dta"
keep if _merge==3
summ
misstable summ
drop _merge
save NHANESdata_PHTHTE_EPHPP_UCPREG_H, replace


***LOAD MERGED DATA FROM SURVEYS D-H 2005-2014
clear all
use NHANESdata_PHTHTE_EPH_PP_UCPREG_D
gen survey="D"
summ 
append using NHANESdata_PHTHTE_EPH_PP_UCPREG_E
replace survey="E" if survey==""
summ
append using NHANESdata_PHTHTE_EPH_PP_UCPREG_F
replace survey="F" if survey==""
summ
append using NHANESdata_PHTHTE_EPH_PP_UCPREG_G
replace survey="G" if survey==""
summ 
append using NHANESdata_PHTHTE_EPHPP_UCPREG_H
replace survey="H" if survey==""
summ
 
*drop exposures not measured in all 5 surveys
drop *mcp* *mnm* *mop* *4to* *opp* *1tb* *3tb* *bpf* *bps* *tlc* *mhnc* *mch*

*drop 4 observations with all exposure measurements missing
drop if urxcop==. & urxcnp==. & urxecp==. & urxmbp==. & urxmc1==. & urxmep==. & urxmhh==. & urxmhp==. & urxmib==. & urxmnp==. & urxmoh==. & urxmzp==. & urxbph==. & urxtrs==. & urxbp3==. & urxppb==. & urxbup==. & urxepb==. & urxmpb==. & urx14d==. & urxdcb==. 
	
*check below LOD values are LOD/sqrt(2) 
local exposures "ecp mep mpb moh mbp mhh bp3 mzp ppb cop mib bph cnp mc1 trs dcb 14d mhp mnp bup epb"
foreach x in `exposures' {
	di "`x'"
	tab urx`x' if urd`x'lc==1
}
	
*rename BPH to BPA, common name
gen urxbpa=urxbph
gen urdbpalc=urdbphlc
drop urxbph urdbphlc

*rename 14D to D14 (as variable names can't start with number)
gen urxd14=urx14d
gen urdd14lc=urd14dlc
drop urx14d urd14dlc

*create creatinine-adjusted exposure vars
local exposures "ecp mep mpb moh mbp mhh bp3 mzp ppb cop mib bpa cnp mc1 trs dcb d14 mhp mnp bup epb"
foreach x in `exposures' {
	gen urx`x'_ucr=urx`x'/urxucr
}

*label new variables
label var urxcop_ucr "Urinary mono(carboxyoctyl) phthalate (ug/g creatinine)"
label var urxcnp_ucr "Urinary mono(carboxynonyl) phthalate (ug/g creatinine)"
label var urxecp_ucr "Urinary mono-2-ethyl-5-carboxypentyl phthalate (ug/g creatinine)"
label var urxmbp_ucr "Urinary mono-n-butyl phthalate (ug/g creatinine)"
label var urxmc1_ucr "Urinary mono-(3-carboxypropyl) phthalate (ug/g creatinine)"
label var urxmep_ucr "Urinary mono-ethyl phthalate (ug/g creatinine)"
label var urxmhh_ucr "Urinary mono-(2-ethyl-5-hydroxyhexyl) phthalate (ug/g creatinine)"
label var urxmhp_ucr "Urinary mono-2-ethyl-5-carboxypentyl phthalate (ug/g creatinine)"
label var urxmib_ucr "Urinary mono-isobutyl phthalate (ug/g creatinine)"
label var urxmnp_ucr "Urinary mono-isononyl phthalate (ug/g creatinine)"
label var urxmoh_ucr "Urinary mono-(2-ethyl-5-oxohexyl) phthalate (ug/g creatinine)"
label var urxmzp_ucr "Urinary mono-benzyl phthalate (ug/g creatinine)"
label var urxbpa_ucr "Urinary bisphenol A (ug/g creatinine)"
label var urxtrs_ucr "Urinary triclosan (ug/g creatinine)" 
label var urxbp3_ucr "Urinary benzophenone-3 (ug/g creatinine)"
label var urxppb_ucr "Urinary propylparaben (ug/g creatinine)"
label var urxbup_ucr "Urinary butylparaben (ug/g creatinine)"
label var urxepb_ucr "Urinary ethylparaben (ug/g creatinine)" 
label var urxmpb_ucr "Urinary methylparaben (ug/g creatinine)" 
label var urxd14_ucr "Urinary 2,5-dichlorophenol (ug/g creatinine)" 
label var urxdcb_ucr "Urinary 2,4-dichlorophenol (ug/g creatinine)" 

*detection rates
summ urdcoplc urdcnplc urdecplc urdmbplc urdmc1lc urdmeplc urdmhhlc urdmhplc urdmiblc urdmnplc urdmohlc urdmzplc urdbpalc urdtrslc urdbp3lc urdppblc urdbuplc urdepblc urdmpblc urdd14lc urddcblc

*keep only exposures with greater than 80% detection rate
local exposures "ecp mep mpb moh mbp mhh bp3 mzp ppb cop mib bpa cnp mc1 trs dcb d14 mhp mnp bup epb"
foreach x in `exposures' {
	summ urd`x'lc, meanonly
	if r(mean)<=0.2 {
		local keeplist "`keeplist' urx`x'_ucr"
	}
}
keep `keeplist'
 
*rename all exposures to 3 letter uppercase abbreviation
qui ds urx*
local varlist `"`r(varlist)'"'
foreach var in `varlist' {	
	local varabb=upper(substr("`var'",4,3))
	rename `var' `varabb'
}


*save
save NHANES_ExposureData_17, replace
export delimited NHANES_ExposureData_17, replace

