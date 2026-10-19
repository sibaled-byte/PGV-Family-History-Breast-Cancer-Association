* Plotting data are aggregate subtype counts obtained from the main analysis.
* This script does not process individual-level records.

* eFigure 3
* TNBC according to BRCA1/2 carrier status and menopausal status

clear

input menopause BRCA_new numerator denominator
0 0 34 366
0 1  5  22
1 0 51 362
1 1  0   4
end

label define menopause_lbl ///
    0 "Premenopausal" ///
    1 "Postmenopausal"

label values menopause menopause_lbl



* CALCULATE PERCENTAGES


* numerator   = number of TNBC cases
* denominator = number of cases with complete subtype information

gen double percentage = 100 numerator / denominator



* BAR POSITIONS


* Premenopausal: 1 and 2
* Postmenopausal: 4 and 5

gen double xpos = 1 + BRCA_new + 3 menopause



* BAR LABELS


* Display percentages as whole numbers
* Actual bar heights retain the exact percentages

gen str8 pct_label = ///
    strtrim(string(percentage, "%5.0f")) + "%"

* Position percentage labels above bars
gen double label_y = percentage + 0.8

* Make the 0% label visible
replace label_y = 0.8 if percentage == 0



* GRAPH


twoway ///
    (bar percentage xpos if BRCA_new == 0, ///
        base(0) ///
        barwidth(.62) ///
        fcolor("190 190 190") ///
        lcolor("190 190 190")) ///
    (bar percentage xpos if BRCA_new == 1, ///
        base(0) ///
        barwidth(.62) ///
        fcolor("166 77 94") ///
        lcolor("166 77 94")) ///
    (scatter label_y xpos, ///
        msymbol(none) ///
        mlabel(pct_label) ///
        mlabposition(0) ///
        mlabsize(small) ///
        mlabcolor(black)), ///
    xlabel( ///
    1 `""{it:BRCA1/2}-" "n = 366""' ///
    2 `""{it:BRCA1/2}+" "n = 22""' ///
    4 `""{it:BRCA1/2}-" "n = 362""' ///
    5 `""{it:BRCA1/2}+" "n = 4""' ///
        , ///
        angle(0) ///
        labsize(small) ///
        noticks nogrid) ///
    ylabel(0(5)30, ///
        angle(horizontal) ///
        labsize(small) ///
        nogrid) ///
    yscale(range(0 30)) ///
    xscale(range(.5 5.5)) ///
    ytitle( ///
        "Triple-negative breast cancer, %", ///
        size(medsmall)) ///
    xtitle( ///
        "Premenopausal                                      Postmenopausal", ///
        size(medsmall) ///
        margin(t=4)) ///
    legend(off) ///
    title("") ///
    graphregion( ///
        color(white) ///
        margin(l=5 r=5 t=3 b=3)) ///
    plotregion( ///
        color(white) ///
        margin(l=3 r=3 t=2 b=0)) ///
    xsize(8) ///
    ysize(5) ///
    name(TNBC_BRCA_menopause, replace)



* EXPORT


graph export "TNBC_BRCA_menopause_corrected.png", ///
    width(2400) replace

graph export "TNBC_BRCA_menopause_corrected.pdf", ///
    replace



* CHECK VALUES


list menopause BRCA_new numerator denominator percentage, noobs
