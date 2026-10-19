version 18
clear all
set more off

 
* eFIGURE 2
* Overall breast cancer-related PGV prevalence by age at diagnosis
* One ordinary bar per age group; no stacked gene-distribution segments
 

use "study_input.dta", clear

 
* 1. DEFINE PGV GROUPS
* Definitions retained from the original do-file.
 

* Any PGV across all 27 sequenced genes
egen any_pgv = rowmax( ///
    APC ATM BARD1 BMPR1A BRCA1 BRCA2 BRIP1 CDH1 CDK4 ///
    CDKN2A CHEK2 EPCAM HOXB13 MLH1 MSH2 MSH6 MUTYH NBN NF1 PALB2 PMS2 ///
    PTEN RAD51C RAD51D SMAD4 STK11 TP53 ///
    ) if sampleqc=="QC clear"

* Breast cancer-related genes
egen bc_pgv = rowmax( ///
    ATM BRCA1 BRCA2 CHEK2 NF1 PALB2 PTEN TP53 ///
    ) if sampleqc=="QC clear"

gen byte BC_new = .
replace BC_new = 1 if bc_pgv==1 & sampleqc=="QC clear"
* Reference category excludes carriers of PGVs in other sequenced genes.
replace BC_new = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BC_new "Breast cancer-related PGV carrier"

* Homologous recombination repair genes
egen hrr_pgv = rowmax( ///
    ATM BARD1 BRCA1 BRCA2 BRIP1 CHEK2 PALB2 RAD51C RAD51D ///
    ) if sampleqc=="QC clear"

gen byte HR_new = .
replace HR_new = 1 if hrr_pgv==1 & sampleqc=="QC clear"
replace HR_new = 0 if any_pgv==0 & sampleqc=="QC clear"
label var HR_new "HRR PGV carrier"

* BRCA1/2
gen byte BRCA_new = .
replace BRCA_new = 1 if (BRCA1==1 | BRCA2==1) & sampleqc=="QC clear"
replace BRCA_new = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BRCA_new "BRCA1/2 PGV carrier"

* HRR excluding BRCA1/2
egen hrr_no_brca_pgv = rowmax( ///
    ATM BARD1 BRIP1 CHEK2 PALB2 RAD51C RAD51D ///
    ) if sampleqc=="QC clear"

gen byte HR_noBRCA = .
replace HR_noBRCA = 1 if hrr_no_brca_pgv==1 & sampleqc=="QC clear"
replace HR_noBRCA = 0 if any_pgv==0 & sampleqc=="QC clear"
label var HR_noBRCA "HRR PGV carrier excluding BRCA1/2"

* BRCA1 separately
gen byte BRCA1_only = .
replace BRCA1_only = 1 if BRCA1==1 & sampleqc=="QC clear"
replace BRCA1_only = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BRCA1_only "BRCA1 PGV carrier"

* BRCA2 separately
gen byte BRCA2_only = .
replace BRCA2_only = 1 if BRCA2==1 & sampleqc=="QC clear"
replace BRCA2_only = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BRCA2_only "BRCA2 PGV carrier"

 
* 2. DEFINE ANALYTIC SAMPLE
* Exclusions retained exactly from the original do-file.
 

drop if sex_check=="NG"
keep if CaCo_status==0 | breast==1
keep if sampleqc=="QC clear"
keep if !missing(breast, breast_family, menopause, fvage, hpc)
drop if menopause==9

display "--"
display "FINAL ANALYTIC SAMPLE"
display "-----"

tab breast, missing
tab breast_family, missing
tab menopause, missing

 
* 3. CREATE AGE GROUP
 

capture drop agegrp

gen byte agegrp = .
replace agegrp = 1 if fvage < 40
replace agegrp = 2 if fvage >=40 & fvage <50
replace agegrp = 3 if fvage >=50 & fvage <60
replace agegrp = 4 if fvage >=60 & fvage <70
replace agegrp = 5 if fvage >=70 & !missing(fvage)

label define agegrp_lbl ///
    1 "<40" ///
    2 "40–49" ///
    3 "50–59" ///
    4 "60–69" ///
    5 "≥70", replace

label values agegrp agegrp_lbl

tab agegrp, missing

 
* 4. CHECK UNIQUE BC-RELATED PGV PREVALENCE
 

display "-----"
display "UNIQUE BC-RELATED PGV PREVALENCE AMONG BREAST CANCER CASES"
display "--"

tab agegrp bc_pgv if breast==1, row
tab agegrp bc_pgv if breast==1

 
* 5. COCHRAN-ARMITAGE TREND TEST
* Calculation retained from the original do-file.
 

preserve

keep if breast==1
keep if inlist(agegrp,1,2,3,4,5)

collapse (sum) x=bc_pgv (count) n=breast, by(agegrp)

gen score = agegrp

summ x
local X = r(sum)

summ n
local Ntotal = r(sum)

local p_pool = `X'/`Ntotal'

summ score [fw=n]
local wbar = r(mean)

gen num_part = score*(x - n*`p_pool')
summ num_part
local numerator = r(sum)

gen var_part = n*(score - `wbar')^2
summ var_part
local varsum = r(sum)

local z = `numerator' / ///
    sqrt(`p_pool'*(1-`p_pool')*`varsum')

local p = 2*normal(-abs(`z'))

* JAMA-style P-value display
if `p' < 0.001 {
    local ptrend "<.001"
}
else {
    local ptrend : display %4.3f `p'
    local ptrend = subinstr("`ptrend'","0.","." ,1)
}

display "Cochran-Armitage P for trend = `ptrend'"

restore

 
* 6. CREATE PLOT DATA
* Overall ("All") plus age-specific prevalence
 

preserve

keep if breast==1
keep if inlist(agegrp,1,2,3,4,5)

* Save age-specific data
tempfile agespecific

* Numerator = BC-related PGV carriers
* Denominator = total breast cancer cases
bysort agegrp: egen n_carrier = total(bc_pgv)
bysort agegrp: gen N = _N

collapse ///
    (mean) bc_pgv ///
    (first) n_carrier N, ///
    by(agegrp)

save `agespecific'

* Create overall ("All") observation
restore
preserve

keep if breast==1
keep if inlist(agegrp,1,2,3,4,5)

count
local Nall = r(N)

quietly count if bc_pgv==1
local nall = r(N)

quietly summarize bc_pgv
local pall = r(mean)*100

* Return to age-specific collapsed data
use `agespecific', clear

* Convert prevalence to %
replace bc_pgv = bc_pgv*100

* Add one observation for "All"
set obs `=_N+1'

replace agegrp = 0 in L
replace bc_pgv = `pall' in L
replace n_carrier = `nall' in L
replace N = `Nall' in L

sort agegrp

label var bc_pgv "Breast cancer-related PGV prevalence (%)"

* Percentage labels above bars

* Position labels above bars
gen label_y = bc_pgv + 0.45


gen str30 pct_label = ///
    string(n_carrier,"%9.0f") + "/" + ///
    string(N,"%9.0f") + " (" + ///
    string(bc_pgv,"%4.1f") + "%)"

list agegrp n_carrier N bc_pgv, noobs

 
* 7. DRAW ORDINARY BAR GRAPH
 

twoway ///
    (bar bc_pgv agegrp, ///
        barwidth(0.62) ///
        color("145 35 55") ///
        lcolor(none)) ///
    (scatter label_y agegrp, ///
        msymbol(none) ///
        mlabel(pct_label) ///
        mlabcolor(black) ///
        mlabsize(vsmall) ///
        mlabposition(0)) ///
    (function y=0, range(-0.4 5.6) lcolor(black) lwidth(vthin)) ///
    , ///
xlabel( ///
    0 "All" ///
    1 "<40" ///
    2 "40–49" ///
    3 "50–59" ///
    4 "60–69" ///
    5 "≥70" ///
    , labsize(small) nogrid) ///
    ylabel(0(2)12, ///
        labsize(small) ///
        angle(0) ///
        nogrid) ///
    yscale(range(0 13) noextend) ///
    xscale(range(-0.4 5.6) noextend) ///
    plotregion(margin(zero) color(white) lcolor(none)) ///
text(11.8 3.5 "{it:P} for trend = `ptrend'", ///
    size(small) place(w)) ///
ytitle("Breast cancer-related PGV prevalence (%)", ///
    size(small)) ///
xtitle("Age at diagnosis, years", ///
    size(small) margin(t=4)) ///
legend(off) ///
graphregion(color(white)) ///
bgcolor(white) ///
xsize(5.8) ///
ysize(3.4)

graph export "breast_PGV_prevalence_by_age_CA.tif", ///
    replace width(3600)

restore

 
* 8. FINAL CHECKS IN FULL ANALYTIC SAMPLE
 

display "----"
display "FINAL CHECKS"
display "----"

tab breast_family if breast==1 & menopause==0 & bc_pgv==1, missing
tab breast_family if breast==1 & menopause==0 & BC_new==1, missing

count if breast==1 & menopause==0 & bc_pgv==1 & breast_family==0
count if breast==1 & menopause==0 & bc_pgv==1 & breast_family==1

* BC-related PGV prevalence by age at diagnosis
tab agegrp bc_pgv if breast==1, row

* Exact counts by age group
tab agegrp bc_pgv if breast==1
