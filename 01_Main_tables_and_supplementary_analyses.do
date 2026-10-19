clear all
set more off

* Joint associations of breast cancer-related PGVs and breast cancer family history
* HERPACC breast cancer case-control analysis
* Outcome: breast cancer (breast)
* Family history: breast cancer in a parent or sibling (breast_family)
* Menopause: 0 = premenopausal, 1 = postmenopausal

capture log close
log using "analysis_log.log", text replace
use "study_input.dta", clear


* Define whether a participant carries any PGV in the 27 sequenced genes.
* The common reference group comprises participants without PGVs in any of the 27 genes.
* Carriers of PGVs outside the gene group under analysis are not classified as noncarriers
* for that group; their group indicator remains missing unless they carry a PGV in the group.
egen any_pgv = rowmax(APC ATM BARD1 BMPR1A BRCA1 BRCA2 BRIP1 CDH1 CDK4 ///
    CDKN2A CHEK2 EPCAM HOXB13 MLH1 MSH2 MSH6 MUTYH NBN NF1 PALB2 PMS2 ///
    PTEN RAD51C RAD51D SMAD4 STK11 TP53) if sampleqc=="QC clear"

* Breast cancer-related genes
egen bc_pgv = rowmax(ATM BRCA1 BRCA2 CHEK2 NF1 PALB2 PTEN TP53) ///
    if sampleqc=="QC clear"

gen byte BC_new = .
replace BC_new = 1 if bc_pgv==1 & sampleqc=="QC clear"
replace BC_new = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BC_new "Breast cancer-related PGV carrier"

* Homologous recombination repair genes
egen hrr_pgv = rowmax(ATM BARD1 BRCA1 BRCA2 BRIP1 CHEK2 PALB2 RAD51C RAD51D) ///
    if sampleqc=="QC clear"

gen byte HR_new = .
replace HR_new = 1 if hrr_pgv==1 & sampleqc=="QC clear"
replace HR_new = 0 if any_pgv==0 & sampleqc=="QC clear"
label var HR_new "HRR PGV carrier"

* BRCA1/2
gen byte BRCA_new = .
replace BRCA_new = 1 if (BRCA1==1 | BRCA2==1) & sampleqc=="QC clear"
replace BRCA_new = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BRCA_new "BRCA1/2 PGV carrier"

* Sensitivity analysis: HRR genes excluding BRCA1 and BRCA2
egen hrr_no_brca_pgv = rowmax(ATM BARD1 BRIP1 CHEK2 PALB2 RAD51C RAD51D) ///
    if sampleqc=="QC clear"

gen byte HR_noBRCA = .
replace HR_noBRCA = 1 if hrr_no_brca_pgv==1 & sampleqc=="QC clear"
replace HR_noBRCA = 0 if any_pgv==0 & sampleqc=="QC clear"
label var HR_noBRCA "HRR PGV carrier excluding BRCA1/2"

* Exploratory analyses of BRCA1 and BRCA2 separately
gen byte BRCA1_only = .
replace BRCA1_only = 1 if BRCA1==1 & sampleqc=="QC clear"
replace BRCA1_only = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BRCA1_only "BRCA1 PGV carrier"

gen byte BRCA2_only = .
replace BRCA2_only = 1 if BRCA2==1 & sampleqc=="QC clear"
replace BRCA2_only = 0 if any_pgv==0 & sampleqc=="QC clear"
label var BRCA2_only "BRCA2 PGV carrier"

* Analysis sample
drop if sex_check=="NG"
keep if CaCo_status==0 | breast==1
keep if sampleqc=="QC clear"
keep if !missing(breast, breast_family, menopause, fvage, hpc)
drop if menopause==9

* Analytic sample distributions
tab breast, missing
tab breast_family, missing
tab menopause, missing

save "analysis_intermediate.dta", replace




* TABLE 1. BASELINE CHARACTERISTICS OF CASES AND CONTROLS
* Age categories
* <40, 40-49, 50-59, 60-69, >=70 years
 

capture drop age_cat
gen age_cat = .

replace age_cat = 1 if fvage < 40
replace age_cat = 2 if fvage >= 40 & fvage < 50
replace age_cat = 3 if fvage >= 50 & fvage < 60
replace age_cat = 4 if fvage >= 60 & fvage < 70
replace age_cat = 5 if fvage >= 70 & fvage < .

label define age_cat_lbl ///
    1 "<40" ///
    2 "40-49" ///
    3 "50-59" ///
    4 "60-69" ///
    5 ">=70", replace

label values age_cat age_cat_lbl

* BMI categories
* <23.0, 23.0-24.9, >=25.0 kg/m2
 

capture drop bmi_cat
gen bmi_cat = .

replace bmi_cat = 1 if bmi < 23
replace bmi_cat = 2 if bmi >= 23 & bmi < 25
replace bmi_cat = 3 if bmi >= 25 & bmi < .

label define bmi_cat_lbl ///
    1 "<23.0" ///
    2 "23.0-24.9" ///
    3 ">=25.0", replace

label values bmi_cat bmi_cat_lbl

* Alcohol consumption categories
* 0, <23, 23-45, >=46 g/day


capture drop alcohol_cat
gen alcohol_cat = .

replace alcohol_cat = 0 if quant_drink_gr == 0
replace alcohol_cat = 1 if quant_drink_gr > 0 & quant_drink_gr < 23
replace alcohol_cat = 2 if quant_drink_gr >= 23 & quant_drink_gr < 46
replace alcohol_cat = 3 if quant_drink_gr >= 46 & quant_drink_gr < .

label define alcohol_cat_lbl ///
    0 "0" ///
    1 "<23" ///
    2 "23-45" ///
    3 ">=46", replace

label values alcohol_cat alcohol_cat_lbl
tab breast alcohol_cat, row chi2

* Age, mean (SD)
tabstat fvage, by(breast) statistics(n mean sd)

* Categorical variables: n (%) and Pearson chi-square P values
tab breast age_cat, row chi2
tab breast menopause, row chi2
tab breast bmi_cat, row chi2
tab breast alcohol_cat, row chi2
tab breast breast_family, row chi2

* carrier status by case-control status
foreach g in BC_new BRCA_new HR_new {
    tab breast `g', chi2
}
* PGV carrier status
foreach g in BC_new BRCA_new HR_new {
    tab breast `g', row chi2
}


* Table 2: Joint associations of PGV carrier status and breast cancer FH
* Adjusted for age and HERPACC recruitment phase
* The # specification gives joint-category ORs relative to FH-negative noncarriers.
* The ## specification includes lower-order terms and is used to assess interaction.
foreach g in BC_new BRCA_new HR_new {
    table breast_family `g' breast, missing
    logit breast i.`g'#i.breast_family c.fvage i.hpc, or
    logit breast i.`g'##i.breast_family c.fvage i.hpc, or
}



* Table 3: Joint associations stratified by menopausal status
* Adjusted for age and HERPACC recruitment phase
* Within each stratum, # reports joint-category ORs and ## fits the factorial interaction model.
foreach g in BC_new BRCA_new HR_new {
    table breast_family menopause breast if `g'==0, missing
    table breast_family menopause breast if `g'==1, missing

    logit breast i.`g'#i.breast_family c.fvage i.hpc if menopause==0, or
    logit breast i.`g'##i.breast_family c.fvage i.hpc if menopause==0, or

    logit breast i.`g'#i.breast_family c.fvage i.hpc if menopause==1, or
    logit breast i.`g'##i.breast_family c.fvage i.hpc if menopause==1, or
}

* Three-way interaction: PGV x FH x menopausal status
foreach g in BC_new BRCA_new HR_new {
    logit breast i.`g'##i.breast_family##i.menopause c.fvage i.hpc, or
    testparm i.`g'#i.breast_family#i.menopause
}


* eTable 2: Participant characteristics according to
* availability of adjustment covariates

capture drop complete_covariates
gen byte complete_covariates = ///
    !missing(fvage, hpc, bmi, quant_drink_gr, py)

label define complete_cov_lbl ///
    0 "At least one missing covariate" ///
    1 "Complete covariates", replace
label values complete_covariates complete_cov_lbl

* Check group sizes
tab complete_covariates, missing


* Continuous characteristics

* Age
tabstat fvage, by(complete_covariates) statistics(n mean sd)
ttest fvage, by(complete_covariates) unequal

* BMI
tabstat bmi, by(complete_covariates) statistics(n mean sd)
ttest bmi, by(complete_covariates) unequal

* Alcohol consumption
tabstat quant_drink_gr, by(complete_covariates) statistics(n median p25 p75)
ranksum quant_drink_gr, by(complete_covariates)

* Smoking exposure
tabstat py, by(complete_covariates) statistics(n median p25 p75)
ranksum py, by(complete_covariates)

* Categorical characteristics
tab complete_covariates breast, row chi2
tab complete_covariates menopause, row chi2
tab complete_covariates breast_family, row chi2
tab complete_covariates hpc, row chi2

* Carrier prevalence using full column denominators
tab complete_covariates bc_pgv, row chi2

gen byte brca12_pgv = (BRCA1 == 1 | BRCA2 == 1) if sampleqc == "QC clear"
tab complete_covariates brca12_pgv, row chi2
drop brca12_pgv

tab complete_covariates hrr_pgv, row chi2



* eTable 3: Unadjusted joint associations by menopausal status

foreach g in BC_new BRCA_new HR_new {

    di "======="
    di "UNADJUSTED ANALYSIS: `g'"
    di "======="

     
    * Premenopausal women
     

    * Cell counts
    table breast_family breast if menopause==0 & `g'==0, missing
    table breast_family breast if menopause==0 & `g'==1, missing

    * Joint ORs:
    * reference = FH-negative noncarriers
    logit breast i.`g'#i.breast_family ///
        if menopause==0, or

    * PGV x FH interaction
    logit breast i.`g'##i.breast_family ///
        if menopause==0, or


     
    * Postmenopausal women
     

    * Cell counts
    table breast_family breast if menopause==1 & `g'==0, missing
    table breast_family breast if menopause==1 & `g'==1, missing

    * Joint ORs:
    * reference = FH-negative noncarriers
    logit breast i.`g'#i.breast_family ///
        if menopause==1, or

    * PGV x FH interaction
    logit breast i.`g'##i.breast_family ///
        if menopause==1, or


     
    * Three-way interaction:
    * PGV x FH x menopausal status
    * Includes all corresponding main effects and 2-way interactions
     

    logit breast i.`g'##i.breast_family##i.menopause, or

    testparm i.`g'#i.breast_family#i.menopause
}


* eTable 4: Fully adjusted models
foreach g in BC_new BRCA_new HR_new {
    logit breast i.`g'#i.breast_family c.fvage i.hpc ///
        c.bmi c.quant_drink_gr c.py, or
    logit breast i.`g'##i.breast_family c.fvage i.hpc ///
        c.bmi c.quant_drink_gr c.py, or
}


* eTable 5: Fully adjusted models by menopausal status
foreach g in BC_new BRCA_new HR_new HR_noBRCA BRCA1_only BRCA2_only {

    di "======= `g' ========"

    table breast_family menopause breast if `g'==0, missing
    table breast_family menopause breast if `g'==1, missing

    * Premenopausal: joint ORs
    logit breast i.`g'#i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==0, or

    * Premenopausal: interaction
    logit breast i.`g'##i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==0, or

    * Postmenopausal: joint ORs
    logit breast i.`g'#i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==1, or

    * Postmenopausal: interaction
    logit breast i.`g'##i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==1, or

    * Three-way interaction
    logit breast i.`g'##i.breast_family##i.menopause ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py, or

    capture noisily testparm i.`g'#i.breast_family#i.menopause
}


* eTable 6: Fully adjusted models stratified by age (<50 and >=50 years)

capture drop age50
gen byte age50 = fvage >= 50 if !missing(fvage)

label define age50_lbl 0 "<50 years" 1 ">=50 years", replace
label values age50 age50_lbl

foreach g in BC_new BRCA_new HR_new HR_noBRCA BRCA1_only BRCA2_only {

    foreach a in 0 1 {

        * Joint ORs
        logit breast i.`g'#i.breast_family ///
            c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
            if age50==`a', or

        * Interaction
        logit breast i.`g'##i.breast_family ///
            c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
            if age50==`a', or
    }
}




* eTable 7: Firth penalized logistic regression sensitivity analyses
foreach g in BC_new BRCA_new HR_new HR_noBRCA BRCA1_only BRCA2_only{

    * Premenopausal: joint ORs for table
    firthlogit breast i.`g'#i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==0, or

    * Premenopausal: interaction
    firthlogit breast i.`g'##i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==0, or

    * Postmenopausal: joint ORs for table
    firthlogit breast i.`g'#i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==1, or

    * Postmenopausal: interaction
    firthlogit breast i.`g'##i.breast_family ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py ///
        if menopause==1, or

    * Three-way interaction
    firthlogit breast i.`g'##i.breast_family##i.menopause ///
        c.fvage i.hpc c.bmi c.quant_drink_gr c.py, or
}



* eTable 8: Breast cancer subtype distribution by PGV carrier status
* and menopausal status

preserve

use "subtype_input.dta", clear

capture drop _merge

local final_indicators ER_p_HER2_n ER_p_HER2_p ER_n_HER2_p ER_n_HER2_n ///
                       ER_p_HER2_unknown ER_n_HER2_unknown

foreach v of local final_indicators {
    capture confirm numeric variable `v'
    if _rc {
        destring `v', replace force
    }
}

gen subtype_record = .
replace subtype_record = 1 if ER_p_HER2_n == 1
replace subtype_record = 2 if ER_p_HER2_p == 1
replace subtype_record = 3 if ER_n_HER2_p == 1
replace subtype_record = 4 if ER_n_HER2_n == 1
replace subtype_record = 5 if ER_p_HER2_unknown == 1
replace subtype_record = 6 if ER_n_HER2_unknown == 1

label define subtype_lbl ///
    1 "ER+/HER2-" ///
    2 "ER+/HER2+" ///
    3 "ER-/HER2+" ///
    4 "Triple-negative / ER-/HER2-" ///
    5 "ER+/HER2 unknown" ///
    6 "ER-/HER2 unknown", replace

label values subtype_record subtype_lbl

gen subtype4_record = subtype_record if inlist(subtype_record, 1, 2, 3, 4)
label values subtype4_record subtype_lbl

* Select earliest available complete subtype for each patient

capture confirm variable 手術日
if !_rc {

    gen complete_subtype_date = 手術日 if subtype4_record < .

    bysort acc: egen first_available_subtype_date = min(complete_subtype_date)

    gen selected_subtype_record = ///
        (手術日 == first_available_subtype_date & subtype4_record < .)

    bysort acc: egen selected_subtype_min = ///
        min(cond(selected_subtype_record == 1, subtype4_record, .))

    bysort acc: egen selected_subtype_max = ///
        max(cond(selected_subtype_record == 1, subtype4_record, .))

    gen subtype_final = selected_subtype_min ///
        if selected_subtype_min == selected_subtype_max
}
else {
    bysort acc: egen subtype_final = min(subtype4_record)
}

label values subtype_final subtype_lbl

* Keep one record per patient

bysort acc: keep if _n == 1

keep acc subtype_final

tempfile subtype_data
save `subtype_data'

restore


* Merge subtype information into the analysis dataset

merge m:1 acc using `subtype_data', keep(master match) nogen

gen subtype = subtype_final
label values subtype subtype_lbl

gen TNBC_final = (subtype == 4) if subtype < .


* eTable 8
* Premenopausal women

display "Premenopausal: Breast cancer-related PGV"
tab subtype BC_new ///
    if breast == 1 & menopause == 0 & subtype < . & !missing(BC_new), ///
    column exact

display "Premenopausal: BRCA1/2"
tab subtype BRCA_new ///
    if breast == 1 & menopause == 0 & subtype < . & !missing(BRCA_new), ///
    column exact


* Postmenopausal women

display "Postmenopausal: Breast cancer-related PGV"
tab subtype BC_new ///
    if breast == 1 & menopause == 1 & subtype < . & !missing(BC_new), ///
    column exact

display "Postmenopausal: BRCA1/2"
tab subtype BRCA_new ///
    if breast == 1 & menopause == 1 & subtype < . & !missing(BRCA_new), ///
    column exact


* TNBC vs all other complete subtypes
* Fisher exact P values reported in eTable 8

display "Premenopausal: TNBC by Breast cancer-related PGV"
tab TNBC_final BC_new ///
    if breast == 1 & menopause == 0 & subtype < . & !missing(BC_new), ///
    column exact

display "Premenopausal: TNBC by BRCA1/2"
tab TNBC_final BRCA_new ///
    if breast == 1 & menopause == 0 & subtype < . & !missing(BRCA_new), ///
    column exact

display "Postmenopausal: TNBC by Breast cancer-related PGV"
tab TNBC_final BC_new ///
    if breast == 1 & menopause == 1 & subtype < . & !missing(BC_new), ///
    column exact

display "Postmenopausal: TNBC by BRCA1/2"
tab TNBC_final BRCA_new ///
    if breast == 1 & menopause == 1 & subtype < . & !missing(BRCA_new), ///
    column exact
* eTable 9: Background characteristics according to subtype-data availability

preserve

keep if breast == 1

count
assert _N == 1282

capture drop subtype_available
gen byte subtype_available = (subtype < .)
label define subtype_avail_lbl 0 "Subtype unavailable" 1 "Subtype available", replace
label values subtype_available subtype_avail_lbl

tab subtype_available, missing

* Continuous characteristics

tabstat fvage, by(subtype_available) statistics(n mean sd median p25 p75)
ttest fvage, by(subtype_available)

tabstat bmi, by(subtype_available) statistics(n mean sd median p25 p75)
ttest bmi, by(subtype_available)

tabstat quant_drink_gr, by(subtype_available) statistics(n mean sd median p25 p75)
ttest quant_drink_gr, by(subtype_available)

tabstat py, by(subtype_available) statistics(n mean sd median p25 p75)
ttest py, by(subtype_available)

* Categorical characteristics

tab subtype_available menopause, row chi2
tab subtype_available breast_family, row chi2
tab subtype_available hpc, row chi2

* Breast cancer-related PGV status

capture drop BCcarrier
gen byte BCcarrier = (bc_pgv == 1) if sampleqc == "QC clear"
label define BCcarrier_lbl 0 "No" 1 "Yes", replace
label values BCcarrier BCcarrier_lbl

tab subtype_available BCcarrier, row chi2

* BRCA1/2 PGV status

capture drop BRCAcarrier
gen byte BRCAcarrier = (BRCA1 == 1 | BRCA2 == 1) if sampleqc == "QC clear"
label define BRCAcarrier_lbl 0 "No" 1 "Yes", replace
label values BRCAcarrier BRCAcarrier_lbl

tab subtype_available BRCAcarrier, row chi2

restore






log close
