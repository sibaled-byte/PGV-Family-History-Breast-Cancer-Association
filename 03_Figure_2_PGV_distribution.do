
* Figure 2
* PGV group distribution among PGV-positive
* side-by-side bars with n/N (%) labels
* premenopausal breast cancer cases by breast cancer family history

clear all
set more off

use "analysis_intermediate.dta", clear

preserve

 
* 1. Restrict to BC-related PGV-positive
*    premenopausal breast cancer cases
 

keep if breast == 1
keep if bc_pgv == 1
keep if menopause == 0
keep if inlist(breast_family, 0, 1)

* Check denominators
tab breast_family

* Expected:
* FH- = 29
* FH+ = 12
* Total = 41


 
* 2. Create mutually exclusive PGV groups
*
* 1 = BRCA1/2
* 2 = HRR excluding BRCA1/2
* 3 = Other breast cancer-related PGVs
 

gen byte BRCA12 = (BRCA1 == 1 | BRCA2 == 1)

gen byte HRRstar = 0
replace HRRstar = 1 if BRCA12 == 0 & ///
    (ATM == 1 | CHEK2 == 1 | PALB2 == 1)

gen byte OtherPGV = 0
replace OtherPGV = 1 if BRCA12 == 0 & HRRstar == 0 & ///
    (NF1 == 1 | PTEN == 1 | TP53 == 1)

* Confirm categories are mutually exclusive and exhaustive
assert BRCA12 + HRRstar + OtherPGV == 1


* 3. Create one categorical PGV-group variable
 

gen byte pgv_group = .

replace pgv_group = 1 if BRCA12 == 1
replace pgv_group = 2 if HRRstar == 1
replace pgv_group = 3 if OtherPGV == 1

label define pgv_lbl ///
    1 "BRCA1/2" ///
    2 "HRR*" ///
    3 "Other PGVs", replace

label values pgv_group pgv_lbl


* 4. Label FH
 

label define fh_lbl ///
    0 "FH-" ///
    1 "FH+", replace

label values breast_family fh_lbl


 
* 5. Confirm raw counts 

tab pgv_group breast_family, column


* 6. Store FH denominators BEFORE collapsing
 

count if breast_family == 0
local N_FHneg = r(N)

count if breast_family == 1
local N_FHpos = r(N)

display "FH- denominator = `N_FHneg'"
display "FH+ denominator = `N_FHpos'"


 
* 7. Collapse to counts by PGV group and FH
 

contract pgv_group breast_family, zero

rename _freq n


 
* 8. Calculate percentages within FH group
 

gen N = .
replace N = `N_FHneg' if breast_family == 0
replace N = `N_FHpos' if breast_family == 1

gen pct = 100 * n / N


 
* 9. Create labels: n/N (%)
*
* Example:
* 23/29 (79%)
* 11/12 (92%)
 

gen str20 barlabel = ///
    string(n, "%2.0f") + "/" + ///
    string(N, "%2.0f") + " (" + ///
    string(pct, "%2.0f") + "%)"


 
* 10. Create positions for side-by-side bars
 

gen xpos = .

* BRCA1/2
replace xpos = 0.85 if pgv_group == 1 & breast_family == 0
replace xpos = 1.15 if pgv_group == 1 & breast_family == 1

* HRR*
replace xpos = 1.85 if pgv_group == 2 & breast_family == 0
replace xpos = 2.15 if pgv_group == 2 & breast_family == 1

* Other PGVs
replace xpos = 2.85 if pgv_group == 3 & breast_family == 0
replace xpos = 3.15 if pgv_group == 3 & breast_family == 1


 
* 11. Check plotting data
 

sort pgv_group breast_family

list pgv_group breast_family n N pct barlabel xpos, ///
    noobs sepby(pgv_group)


 
* 12. Plot clean side-by-side bars
 

twoway ///
    (bar pct xpos if breast_family == 0 & pct > 0, ///
        barwidth(0.24) ///
        base(0) ///
color("39 91 135") ///
lcolor("39 91 135")) ///
    (bar pct xpos if breast_family == 1 & pct > 0, ///
        barwidth(0.24) ///
        base(0) ///
color("166 30 68") ///
lcolor("166 30 68")) ///
    (scatter pct xpos if breast_family == 0 & pct > 0, ///
        msymbol(none) ///
        mlabel(barlabel) ///
        mlabcolor(black) ///
        mlabsize(vsmall) ///
        mlabposition(12)) ///
    (scatter pct xpos if breast_family == 1 & pct > 0, ///
        msymbol(none) ///
        mlabel(barlabel) ///
        mlabcolor(black) ///
        mlabsize(vsmall) ///
        mlabposition(12)) ///
    (scatter pct xpos if breast_family == 1 & pct == 0, ///
        msymbol(none) ///
        mlabel(barlabel) ///
        mlabcolor(black) ///
        mlabsize(vsmall) ///
        mlabposition(12)) ///
    , ///
xlabel( ///
    1 "{it:BRCA1/2}" ///
    2 "HRR*" ///
    3 "Other PGVs", ///
    labsize(small) ///
    noticks ///
    nogrid) ///
    ylabel(0 "0%" 25 "25%" 50 "50%" 75 "75%" 100 "100%", ///
        angle(horizontal) labsize(small) nogrid) ///
    yscale(range(0 105)) ///
    xscale(range(0.55 3.45)) ///
    yline(0, lcolor(black) lwidth(vthin)) ///
    ytitle("Cases, %", size(small)) ///
    xtitle("") ///
    legend(order(1 "FH-" 2 "FH+") ///
        rows(1) ///
        position(6) ///
        size(small) ///
        region(lcolor(none))) ///
    plotregion( ///
        margin(l=2 r=2 b=0 t=3) ///
        color(white) ///
        lcolor(none)) ///
    graphregion( ///
        margin(l=3 r=3 b=2 t=2) ///
        color(white)) ///
    bgcolor(white) ///
    xsize(5.5) ///
    ysize(3.6) ///
    name(Figure2, replace)

 
* 13. Export
 

graph export ///
    "Figure2_PGV_distribution_premenopausal_JAMA.tif", ///
    replace width(3000)

graph export ///
    "Figure2_PGV_distribution_premenopausal_JAMA.png", ///
    replace width(3000)


restore
