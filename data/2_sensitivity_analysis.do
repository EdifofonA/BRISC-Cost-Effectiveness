

*** ONE-WAY SENSITIVITY ANALYSIS ***************************
use Analysis_imputed_data, clear
keep if Rnd == 1
keep treat cost cost_syrp_low cost_syrp_high cost_mnp_low cost_mnp_high cost_prgm_low cost_prgm_high cost_visit_low cost_visit_high cost_bed_low cost_bed_high cost_tx_low cost_tx_high cost_hs* daly* _mi*

foreach var of varlist cost* daly* {
	mi passive: by treat: egen mean_`var' = mean(`var')
}

collapse mean_cost* mean_daly*, by(treat) 
xpose, clear varname
drop in 1
rename v1 dn
rename v2 mnp
rename v3 syrp
order _varname
gen mnp_dn  = mnp  - dn
gen syrp_dn = syrp - dn
drop syrp mnp syrp


generate icer_mnp_base = mnp_dn[1] / mnp_dn[27] in 1/13
replace icer_mnp_base  = mnp_dn[14] / mnp_dn[27] in 14/26

generate icer_syrp_base = syrp_dn[1] / syrp_dn[27] in 1/13
replace icer_syrp_base  = syrp_dn[14] / syrp_dn[27] in 14/26

generate icer_mnp_low   = .
generate icer_mnp_high  = .
generate icer_syrp_low  = .
generate icer_syrp_high = .

replace icer_mnp_low  = (mnp_dn[4]/mnp_dn[27])  in 3
replace icer_mnp_low  = (mnp_dn[6]/mnp_dn[27])  in 4
replace icer_mnp_low  = (mnp_dn[8]/mnp_dn[27])  in 1
replace icer_mnp_low  = (mnp_dn[10]/mnp_dn[27]) in 2
replace icer_mnp_low  = (mnp_dn[1]/mnp_dn[28])  in 5
replace icer_mnp_high = (mnp_dn[5]/mnp_dn[27])  in 3
replace icer_mnp_high = (mnp_dn[7]/mnp_dn[27])  in 4
replace icer_mnp_high = (mnp_dn[9]/mnp_dn[27])  in 1
replace icer_mnp_high = (mnp_dn[11]/mnp_dn[27]) in 2
replace icer_mnp_high = (mnp_dn[1]/mnp_dn[29])  in 5

replace icer_syrp_low  = (syrp_dn[2]/syrp_dn[27])  in 3
replace icer_syrp_low  = (syrp_dn[6]/syrp_dn[27])  in 4
replace icer_syrp_low  = (syrp_dn[8]/syrp_dn[27])  in 1
replace icer_syrp_low  = (syrp_dn[10]/syrp_dn[27]) in 2
replace icer_syrp_low  = (syrp_dn[12]/syrp_dn[28]) in 5
replace icer_syrp_high = (syrp_dn[3]/syrp_dn[27])  in 3
replace icer_syrp_high = (syrp_dn[7]/syrp_dn[27])  in 4
replace icer_syrp_high = (syrp_dn[9]/syrp_dn[27])  in 1
replace icer_syrp_high = (syrp_dn[11]/syrp_dn[27]) in 2
replace icer_syrp_high = (syrp_dn[1]/syrp_dn[29])  in 5

replace icer_mnp_low  = (mnp_dn[17]/mnp_dn[27]) in 8
replace icer_mnp_low  = (mnp_dn[19]/mnp_dn[27]) in 9
replace icer_mnp_low  = (mnp_dn[21]/mnp_dn[27]) in 6
replace icer_mnp_low  = (mnp_dn[23]/mnp_dn[27]) in 7
replace icer_mnp_low  = (mnp_dn[14]/mnp_dn[28]) in 10
replace icer_mnp_high = (mnp_dn[18]/mnp_dn[27]) in 8
replace icer_mnp_high = (mnp_dn[20]/mnp_dn[27]) in 9
replace icer_mnp_high = (mnp_dn[22]/mnp_dn[27]) in 6
replace icer_mnp_high = (mnp_dn[24]/mnp_dn[27]) in 7
replace icer_mnp_high = (mnp_dn[14]/mnp_dn[29]) in 10

replace icer_syrp_low  = (syrp_dn[15]/syrp_dn[27]) in 8
replace icer_syrp_low  = (syrp_dn[19]/syrp_dn[27]) in 9
replace icer_syrp_low  = (syrp_dn[21]/syrp_dn[27]) in 6
replace icer_syrp_low  = (syrp_dn[23]/syrp_dn[27]) in 7
replace icer_syrp_low  = (syrp_dn[12]/syrp_dn[28]) in 10
replace icer_syrp_high = (syrp_dn[16]/syrp_dn[27]) in 8
replace icer_syrp_high = (syrp_dn[20]/syrp_dn[27]) in 9
replace icer_syrp_high = (syrp_dn[22]/syrp_dn[27]) in 6
replace icer_syrp_high = (syrp_dn[24]/syrp_dn[27]) in 7
replace icer_syrp_high = (syrp_dn[14]/syrp_dn[29]) in 10


gen input = _n
order input icer_mnp_base icer_mnp_low icer_mnp_high
replace icer_mnp_base  = icer_mnp_base[14] in 6/10
replace icer_syrp_base = icer_syrp_base[14] in 6/10
keep in 1/10
replace input = input - 5 in 6/10
drop _varname dn mnp_dn syrp_dn


label define INPUT 3 "Cost of active iron agent ({&plusminus} 50%)" 4 "Program delivery cost ({&plusminus} 50%)" 1 "Cost per outpatient visit (95% UI)" 2 "Cost per inpatient bed day (95% UI)" 5 "Anemia disability weight (95% UI)" 
label val input INPUT

generate perspective = 1 in 1/5
replace perspective  = 2 in 6/10
label define PERSPECTIVE 1 "Health system" 2 "Societal"
label values perspective PERSPECTIVE

save Analysis_imputed_data_owsa, replace
set scheme s1color
set cformat %9.1f

use Analysis_imputed_data_owsa, clear
set scheme s1color

** health system perspective: mnp vs placebo
summarize icer_mnp_base in 1/5
local icer = round(r(mean))
twoway (rbar icer_mnp_base icer_mnp_low input, horizontal  bcolor(dkorange%80) barwidth(0.4)) ///
(rbar icer_mnp_base icer_mnp_high input, horizontal bcolor(emerald) barwidth(0.4) xline(`r(mean)', lcolor(gs0) lwidth(thin))) if perspective == 1, ///
subtitle("A. MNPs: ICER  = $ `icer'", size(medium)) ///
xsize(16) ysize(8) ///
ylabel(,           labsize(medsmall) nogrid noticks glcolor(black%3) val angle(horizontal))  ///
xlabel(0(1000)4000, labsize(small) grid glcolor(black%3)) ///
xtitle(, margin(none) size(small)) ///
ytitle(" ", margin(none) ) ///
graphregion(lcolor(white) margin(medium)) ///
plotregion(lwidth(none) fcolor(white)) ///
legend(ring(0) label(1 "Low value") label(2 "High value") rows(2) pos(4) size(small) symxsize(2) symysize(2) region(lwidth(none) color(white%0) )) ///
yscale(noline titlegap(1) lwidth(none)) xscale(titlegap(1) lwidth(medthin)) ///
saving(owsa_mnp_hs, replace)

** health system perspective: mnp vs placebo
summarize icer_syrp_base in 1/5
local icer = round(r(mean))
twoway (rbar icer_syrp_base icer_syrp_low input,  horizontal  bcolor(sandb%80) barwidth(0.4)) ///
(rbar icer_syrp_base icer_syrp_high input, horizontal bcolor(emerald) barwidth(0.4) xline(`r(mean)', lcolor(gs0) lwidth(thin))) if perspective == 1, ///
subtitle(" " "B. Iron supplements: ICER  = $ `icer'", size(medium)) ///
xsize(16) ysize(8) ///
ylabel(,           labsize(medsmall) nogrid noticks glcolor(black%3) val angle(horizontal))  ///
xlabel(0(1000)4000, labsize(small) grid glcolor(black%3)) ///
xtitle(, margin(none) size(small)) ///
ytitle(" ", margin(none) ) ///
graphregion(lcolor(white) margin(medium)) ///
plotregion(lwidth(none) fcolor(white)) ///
legend(ring(0) label(1 "Low value") label(2 "High value") rows(2) pos(4) size(small) symxsize(2) symysize(2) region(lwidth(none) color(white%0) )) ///
yscale(noline titlegap(1) lwidth(none)) xscale(titlegap(1) lwidth(medthin)) ///
saving(owsa_syrp_hs, replace)

** societal perspective: mnp vs placebo
summarize icer_mnp_base in 6/10
local icer = round(r(mean))
twoway (rbar icer_mnp_base icer_mnp_low input,  horizontal  bcolor(dkorange%80) barwidth(0.4)) ///
(rbar icer_mnp_base icer_mnp_high input, horizontal bcolor(emerald) barwidth(0.4) xline(`r(mean)', lcolor(gs0) lwidth(thin))) if perspective == 2, ///
subtitle(" " "C. MNPs (societal): ICER  = $ `icer'", size(medium)) ///
xsize(16) ysize(8) ///
ylabel(,           labsize(medsmall) nogrid noticks glcolor(black%3) val angle(horizontal))  ///
xlabel(0(1000)4000, labsize(small) grid glcolor(black%3)) ///
xtitle(, margin(none) size(small)) ///
ytitle(" ", margin(none) ) ///
graphregion(lcolor(white) margin(medium)) ///
plotregion(lwidth(none) fcolor(white)) ///
legend(ring(0) label(1 "Low value") label(2 "High value") rows(2) pos(4) size(small) symxsize(2) symysize(2) region(lwidth(none) color(white%0) )) ///
yscale(noline titlegap(1) lwidth(none)) xscale(titlegap(1) lwidth(medthin)) ///
saving(owsa_mnp, replace)

** societal perspective: mnp vs placebo
summarize icer_syrp_base in 6/10
local icer = round(r(mean))
twoway (rbar icer_syrp_base icer_syrp_low input,  horizontal  bcolor(sandb%80) barwidth(0.4)) ///
(rbar icer_syrp_base icer_syrp_high input, horizontal bcolor(emerald) barwidth(0.4) xline(`r(mean)', lcolor(gs0) lwidth(thin))) if perspective == 2, ///
subtitle(" " "D. Iron supplements (societal): ICER  = $ `icer'", size(medium)) ///
xsize(16) ysize(8) ///
ylabel(,           labsize(medsmall) nogrid noticks glcolor(black%3) val angle(horizontal))  ///
xlabel(0(1000)4000, labsize(small) grid glcolor(black%3)) ///
xtitle("Cost per DALY averted", margin(small) size(medsmall)) ///
ytitle(" ", margin(none) ) ///
graphregion(lcolor(white) margin(medium)) ///
plotregion(lwidth(none) fcolor(white)) ///
legend(ring(0) label(1 "Low value") label(2 "High value") rows(2) pos(4) size(small) symxsize(2) symysize(2) region(lwidth(none) color(white%0) )) ///
yscale(noline titlegap(1) lwidth(none)) xscale(titlegap(1) lwidth(medthin)) ///
saving(owsa_syrp, replace)
 
/*
graph combine "owsa_mnp_hs" "owsa_syrp_hs"  "owsa_mnp"  "owsa_syrp", col(1) xsize(6) ysize(9) scheme(s1color) imargin(2 2 2 2) saving(owsa_combined, replace) graphregion(color(gs15%35))
*/
graph combine "owsa_mnp_hs" "owsa_syrp_hs", col(1) xsize(6) ysize(5) scheme(s1color) imargin(2 2 2 2) saving(owsa_combined, replace) graphregion(color(white))
graph export Graph_imputed_owsa_hs.tif, replace
graph export Graph_imputed_owsa_hs.eps, replace

cap erase owsa_mnp_hs.gph
cap erase owsa_mnp.gph
cap erase owsa_syrp_hs.gph
cap erase owsa_syrp.gph
cap erase owsa_combined.gph



*** TWO-WAY SENSITIVITY ANALYSIS ***************************
** range of values for mnp strategy at 985
postfile contour cost_prgm cost_intv pce using "twsa_mnp_985.dta", replace 
forvalues cost_prgm = 0(0.1)4.1 {
	forvalues cost_intv = 0(0.1)2.6 {
		use Analysis_imputed_data.dta if treat != 2, clear
		mi set M -= 35
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_mnp  = `cost_intv' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		mi passive: gen nmb = daly * 985 - cost_hs
		mi estimate: qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for iron strategy at 985
postfile contour cost_prgm cost_intv pce using "twsa_syrp_985.dta", replace 
forvalues cost_prgm = 0(0.1)4.1 {
	forvalues cost_intv = 0(0.1)2.6 {
		use Analysis_imputed_data.dta if treat != 1, clear
		mi set M -= 35
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_syrp  = `cost_intv' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		mi passive: gen nmb = daly * 985 - cost_hs
		mi estimate: qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for mnp strategy at 200
postfile contour cost_prgm cost_intv pce using "twsa_mnp_200.dta", replace 
forvalues cost_prgm = 0(0.1)4.1  {
	forvalues cost_intv = 0(0.1)2.6 {
		use Analysis_imputed_data.dta if treat != 2, clear
		mi set M -= 35
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_mnp  = `cost_intv' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		mi passive: gen nmb = daly * 200 - cost_hs
		mi estimate: qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for iron strategy at 200
postfile contour cost_prgm cost_intv pce using "twsa_syrp_200.dta", replace 
forvalues cost_prgm = 0(0.1)4.1  {
	forvalues cost_intv = 0(0.1)2.6 {
		use Analysis_imputed_data.dta if treat != 1, clear
		mi set M -= 35
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_syrp  = `cost_intv' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		mi passive: gen nmb = daly * 200 - cost_hs
		mi estimate: qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour



** range of values for mnp strategy at 1970
postfile contour cost_prgm cost_intv pce using "twsa_mnp_1970.dta", replace 
forvalues cost_prgm = 0(0.1)4.1  {
	forvalues cost_intv = 0(0.1)2.6 {
		use Analysis_imputed_data.dta if treat != 2, clear
		mi set M -= 35
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_mnp  = `cost_intv' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		mi passive: gen nmb = daly * 1970 - cost_hs
		mi estimate: qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for iron strategy at 1970
postfile contour cost_prgm cost_intv pce using "twsa_syrp_1970.dta", replace 
forvalues cost_prgm = 0(0.1)4.1  {
	forvalues cost_intv = 0(0.1)2.6 {
		use Analysis_imputed_data.dta if treat != 1, clear
		mi set M -= 35
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_syrp  = `cost_intv' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		mi passive: gen nmb = daly * 1970 - cost_hs
		mi estimate: qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour



use twsa_mnp_200.dta, clear
generate comparison = 1
append using twsa_syrp_200.dta
replace comparison = 2 if comparison ==.
generate lambda = 200
append using twsa_mnp_985.dta
replace comparison = 1 if comparison ==.
append using twsa_syrp_985.dta
replace comparison = 2 if comparison ==.
replace lambda = 985 if lambda == .
append using twsa_mnp_1970.dta
replace comparison = 1 if comparison ==.
append using twsa_syrp_1970.dta
replace comparison = 2 if comparison ==.
replace lambda = 1970 if lambda == .



label define COMPARISON 1 "MNP vs DN" 2 "Iron vs DN"
label values comparison COMPARISON
label variable comparison "Comparison"
label variable pce "Probability cost-effective"
label variable cost_intv "Cost of intervention material"
label variable cost_prgm "Cost of program management"
label variable lambda "Threshold level"
save Analysis_imputed_data_twsa.dta, replace

cap erase twsa_mnp_200.dta
cap erase twsa_syrp_200.dta
cap erase twsa_mnp_985.dta
cap erase twsa_syrp_985.dta


use Analysis_imputed_data_twsa.dta, clear

* Coloured
** Contour plot for MNP vs DN at $200
twoway (contour pce cost_prgm cost_intv, levels(400) crule(linear) scolor(eggshell) ecolor(olive_teal) zlabel(#6)) (pci 0 1.6 4 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 200, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle(" ", margin(zero) size(medlarge)) ///
ytitle("Program delivery cost ($)", margin(small) size(medlarge)) ///
ztitle("", margin(zero) size(medlarge)) ///
subtitle(" " "MNPs, {&lambda}=200", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 1.6 "$ 1.60", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_mnp_200, replace)

** Contour plot for Iron vs DN at $200
twoway (contour pce cost_prgm cost_intv, levels(400) crule(linear) scolor(eggshell) ecolor(sandb) zlabel(#6)) (pci 0 0.6 4 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 200, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle(" ", margin(zero) size(medlarge)) ///
ytitle("", margin(small) size(medlarge)) ///
ztitle(, margin(small) size(medlarge)) ///
subtitle(" " "Iron supplements, {&lambda}=200", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 0.6 "$ 0.63", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_syrp_200, replace)

** Contour plot for MNP vs DN at $985
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(olive_teal) zlabel(#6)) (pci 0 1.6 4 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 985, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle("Cost of MNPs ($)", margin(small) size(medlarge)) ///
ytitle("Program delivery cost ($)", margin(small) size(medlarge)) ///
ztitle("", margin(zero) size(medlarge)) ///
subtitle(" " "MNPs, {&lambda}=985", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 1.6 "$ 1.60", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_mnp_985, replace)

** Contour plot for Iron vs DN at $985
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(sandb) zlabel(#6)) (pci 0 0.6 4 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 985, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle("Cost of iron syrup ($)", margin(small) size(medlarge)) ///
ytitle("", margin(small) size(medlarge)) ///
ztitle(, margin(small) size(medlarge)) ///
subtitle(" " "Iron supplements, {&lambda}=985", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 0.6 "$ 0.63", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_syrp_985, replace)


** Contour plot for MNP vs DN at $1970
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(olive_teal) zlabel(#6)) (pci 0 1.6 4 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 1970, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle("Cost of MNPs ($)", margin(small) size(medium)) ///
ytitle("Program delivery cost ($)", margin(small) size(medium)) ///
ztitle("", margin(zero) size(medium)) ///
subtitle(" " "MNPs: {&lambda}=1970", size(large)) ///
xlabel(0(0.5)2.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
clegend(title("Probability", size(medsmall)) ) ///
text(3.8 1.6 "$ 1.60", color(gs5) placement(e) box size(.5cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_mnp_1970, replace)

** Contour plot for Iron vs DN at $1970
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(sandb) zlabel(#6)) (pci 0 0.6 4 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 1970, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle("Cost of iron syrup ($)", margin(small) size(medium)) ///
ytitle("", margin(small) size(medium)) ///
ztitle("", margin(small) size(medium)) ///
subtitle(" " "Iron supplements, {&lambda}=1970", size(large)) ///
xlabel(0(0.5)2.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
clegend(title("Probability", size(medsmall)) ) ///
text(3.8 0.6 "$ 0.63", color(gs5) placement(e) box size(.5cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_syrp_1970, replace)



graph combine "twsa_mnp_200" "twsa_syrp_200", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined_200, replace) graphregion(color(white))
graph combine "twsa_mnp_985" "twsa_syrp_985", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined_985, replace) graphregion(color(white))
graph combine "twsa_mnp_1970" "twsa_syrp_1970", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined_1970, replace) graphregion(color(white))
graph export Graph_imputed_twsa_1970.tif, replace

graph combine "twsa_combined_200" "twsa_combined_985", col(1) xsize(15) ysize(12.5) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined, replace) plotregion(margin(none))graphregion(color(white))

*graph combine "twsa_combined_200" "twsa_combined_985" "twsa_combined_1970", col(1) xsize(7.5) ysize(9.375) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined2, replace) plotregion(margin(none))graphregion(color(white))


graph export figure2_color.png, replace
graph export figure2_color.eps, replace




cap erase twsa_mnp_200.gph
cap erase twsa_syrp_200.gph
cap erase twsa_mnp_985.gph
cap erase twsa_syrp_985.gph
cap erase twsa_combined_200.gph
cap erase twsa_combined_985.gph
cap erase twsa_mnp_1970.gph
cap erase twsa_syrp_1970.gph
cap erase twsa_combined.gph
cap erase twsa_combined_1970.gph






use Analysis_imputed_data_twsa.dta, clear
* Black and white
** Contour plot for MNP vs DN at $200
twoway (contour pce cost_prgm cost_intv, levels(400) crule(linear) scolor(gs16) ecolor(gs10) zlabel(#6)) (pci 0 1.6 4 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 200, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle(" ", margin(zero) size(medlarge)) ///
ytitle("Program delivery cost ($)", margin(small) size(medlarge)) ///
ztitle("", margin(zero) size(medlarge)) ///
subtitle(" " "MNPs, {&lambda}=200", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 1.6 "$ 1.60", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_mnp_200, replace)

** Contour plot for Iron vs DN at $200
twoway (contour pce cost_prgm cost_intv, levels(400) crule(linear) scolor(gs16) ecolor(gs8) zlabel(#6)) (pci 0 0.6 4 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 200, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle(" ", margin(zero) size(medlarge)) ///
ytitle("", margin(small) size(medlarge)) ///
ztitle(, margin(small) size(medlarge)) ///
subtitle(" " "Iron supplements, {&lambda}=200", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 0.6 "$ 0.63", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_syrp_200, replace)

** Contour plot for MNP vs DN at $985
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(gs16) ecolor(gs10) zlabel(#6)) (pci 0 1.6 4 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 985, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle("Cost of MNPs ($)", margin(small) size(medlarge)) ///
ytitle("Program delivery cost ($)", margin(small) size(medlarge)) ///
ztitle("", margin(zero) size(medlarge)) ///
subtitle(" " "MNPs, {&lambda}=985", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 1.6 "$ 1.60", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_mnp_985, replace)

** Contour plot for Iron vs DN at $985
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(gs16) ecolor(gs10) zlabel(#6)) (pci 0 0.6 4 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 985, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle("Cost of iron syrup ($)", margin(small) size(medlarge)) ///
ytitle("", margin(small) size(medlarge)) ///
ztitle(, margin(small) size(medlarge)) ///
subtitle(" " "Iron supplements, {&lambda}=985", size(large)) ///
xlabel(0(0.5)2.5, labsize(medlarge) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medlarge) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medlarge) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(3.8 0.6 "$ 0.63", color(gs5) placement(e) box size(1.25cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_syrp_985, replace)

graph combine "twsa_mnp_200" "twsa_syrp_200", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined_200, replace) graphregion(color(white))
graph combine "twsa_mnp_985" "twsa_syrp_985", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined_985, replace) graphregion(color(white))
graph combine "twsa_combined_200" "twsa_combined_985", col(1) xsize(15) ysize(12.5) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined, replace) plotregion(margin(none))graphregion(color(white))

graph export figure2_blackwhite.png, replace
graph export figure2_blackwhite.eps, replace


cap erase twsa_mnp_200.gph
cap erase twsa_syrp_200.gph
cap erase twsa_mnp_985.gph
cap erase twsa_syrp_985.gph
cap erase twsa_combined_200.gph
cap erase twsa_combined_985.gph
cap erase twsa_mnp_1970.gph
cap erase twsa_syrp_1970.gph
cap erase twsa_combined.gph
cap erase twsa_combined_1970.gph



*** TWO-WAY SENSITIVITY ANALYSIS FOR EXTENDED CASE  ***************************
** range of values for mnp strategy at 985
postfile contour cost_prgm cost_intv pce using "twsa_ext_mnp_985.dta", replace 
forvalues cost_prgm = 0(0.1)6.1 {
	forvalues cost_intv = 0(0.1)3.6 {
		use Analysis_complete_data.dta if available == 7 & treat != 2, clear
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_mnp  = `cost_intv' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 985 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for iron strategy at 985
postfile contour cost_prgm cost_intv pce using "twsa_ext_syrp_985.dta", replace 
forvalues cost_prgm = 0(0.1)6.1 {
	forvalues cost_intv = 0(0.1)3.6 {
		use Analysis_complete_data.dta if available == 7  & treat != 1, clear
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_syrp  = `cost_intv' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 985 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for mnp strategy at 200
postfile contour cost_prgm cost_intv pce using "twsa_ext_mnp_200.dta", replace 
forvalues cost_prgm = 0(0.1)6.1  {
	forvalues cost_intv = 0(0.1)3.6 {
		use Analysis_complete_data.dta if available == 7 & treat != 2, clear
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_mnp  = `cost_intv' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 200 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for iron strategy at 200
postfile contour cost_prgm cost_intv pce using "twsa_ext_syrp_200.dta", replace 
forvalues cost_prgm = 0(0.1)6.1  {
	forvalues cost_intv = 0(0.1)3.6 {
		use Analysis_complete_data.dta if available == 7  & treat != 1, clear
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_syrp  = `cost_intv' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 200 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour



** range of values for mnp strategy at 1970
postfile contour cost_prgm cost_intv pce using "twsa_ext_mnp_1970.dta", replace 
forvalues cost_prgm = 0(0.1)6.1  {
	forvalues cost_intv = 0(0.1)3.6 {
		use Analysis_complete_data.dta if available == 7 & treat != 2, clear
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_mnp  = `cost_intv' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 1970 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour

** range of values for iron strategy at 1970
postfile contour cost_prgm cost_intv pce using "twsa_ext_syrp_1970.dta", replace 
forvalues cost_prgm = 0(0.1)6.1  {
	forvalues cost_intv = 0(0.1)3.6 {
		use Analysis_complete_data.dta if available == 7  & treat != 1, clear
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_syrp  = `cost_intv' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 1970 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`cost_intv') (`pce')
	}
}
postclose contour



use twsa_ext_mnp_200.dta, clear
generate comparison = 1
append using twsa_ext_syrp_200.dta
replace comparison = 2 if comparison ==.
generate lambda = 200
append using twsa_ext_mnp_985.dta
replace comparison = 1 if comparison ==.
append using twsa_ext_syrp_985.dta
replace comparison = 2 if comparison ==.
replace lambda = 985 if lambda == .
append using twsa_ext_mnp_1970.dta
replace comparison = 1 if comparison ==.
append using twsa_ext_syrp_1970.dta
replace comparison = 2 if comparison ==.
replace lambda = 1970 if lambda == .



label define COMPARISON 1 "MNP vs DN" 2 "Iron vs DN"
label values comparison COMPARISON
label variable comparison "Comparison"
label variable pce "Probability cost-effective"
label variable cost_intv "Cost of intervention material"
label variable cost_prgm "Cost of program management"
label variable lambda "Threshold level"
save Analysis_extended_data_twsa.dta, replace

cap erase twsa_ext_mnp_200.dta
cap erase twsa_ext_syrp_200.dta
cap erase twsa_ext_mnp_985.dta
cap erase twsa_ext_syrp_985.dta
cap erase twsa_ext_mnp_1970.dta
cap erase twsa_ext_syrp_1970.dta


use Analysis_extended_data_twsa.dta, clear

** Contour plot for MNP vs DN at $200
twoway (contour pce cost_prgm cost_intv, levels(400) crule(linear) scolor(eggshell) ecolor(olive_teal) zlabel(#6)) (pci 0 1.6 6.1 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 200, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle(" ", margin(zero) size(medium)) ///
ytitle("Program delivery cost ($)", margin(small) size(medium)) ///
ztitle("", margin(zero) size(medium)) ///
subtitle(" " "MNPs, {&lambda}=200", size(large)) ///
xlabel(0(0.5)3.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)6, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(5.8 1.6 "$ 1.60", color(gs5) placement(e) box size(.4cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_ext_mnp_200, replace)

** Contour plot for Iron vs DN at $200
twoway (contour pce cost_prgm cost_intv, levels(400) crule(linear) scolor(eggshell) ecolor(sandb) zlabel(#6)) (pci 0 0.6 6.1 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 200, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle(" ", margin(zero) size(medium)) ///
ytitle("", margin(small) size(medium)) ///
ztitle(, margin(small) size(medium)) ///
subtitle(" " "Iron supplements, {&lambda}=200", size(large)) ///
xlabel(0(0.5)3.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)6, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(5.8 0.6 "$ 0.63", color(gs5) placement(e) box size(.4cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_ext_syrp_200, replace)

** Contour plot for MNP vs DN at $985
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(olive_teal) zlabel(#6)) (pci 0 1.6 6.1 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 985, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle("Cost of MNPs ($)", margin(small) size(medium)) ///
ytitle("Program delivery cost ($)", margin(small) size(medium)) ///
ztitle("", margin(zero) size(medium)) ///
subtitle(" " "MNPs, {&lambda}=985", size(large)) ///
xlabel(0(0.5)3.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)6, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(5.8 1.6 "$ 1.60", color(gs5) placement(e) box size(.4cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_ext_mnp_985, replace)

** Contour plot for Iron vs DN at $985
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(sandb) zlabel(#6)) (pci 0 0.6 6.1 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 985, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle("Cost of iron syrup ($)", margin(small) size(medium)) ///
ytitle("", margin(small) size(medium)) ///
ztitle(, margin(small) size(medium)) ///
subtitle(" " "Iron supplements, {&lambda}=985", size(large)) ///
xlabel(0(0.5)3.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)6, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
text(5.8 0.6 "$ 0.63", color(gs5) placement(e) box size(.4cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_ext_syrp_985, replace)


** Contour plot for MNP vs DN at $1970
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(olive_teal) zlabel(#6)) (pci 0 1.6 6.1 1.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 1 & lambda == 1970, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medthin)) ///
xtitle("Cost of MNPs ($)", margin(small) size(medium)) ///
ytitle("Program delivery cost ($)", margin(small) size(medium)) ///
ztitle("", margin(zero) size(medium)) ///
subtitle(" " "MNPs: {&lambda}=1970", size(large)) ///
xlabel(0(0.5)2.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
clegend(title("Probability", size(medsmall)) ) ///
text(5.8 1.6 "$ 1.60", color(gs5) placement(e) box size(.4cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_ext_mnp_1970, replace)

** Contour plot for Iron vs DN at $1970
twoway (contour pce cost_prgm cost_intv, levels(300) crule(linear) scolor(eggshell) ecolor(sandb) zlabel(#6)) (pci 0 0.6 6.1 0.6, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) if comparison == 2 & lambda == 1970, ///
graphregion(lcolor(white) fcolor(white)  margin(none)) ///
plotregion(margin(tiny) fcolor(white) lcolor(gs0) lwidth(medium)) ///
xtitle("Cost of iron syrup ($)", margin(small) size(medium)) ///
ytitle("", margin(small) size(medium)) ///
ztitle("", margin(small) size(medium)) ///
subtitle(" " "Iron supplements, {&lambda}=1970", size(large)) ///
xlabel(0(0.5)2.5, labsize(medsmall) format(%4.1f) nogrid)  ///
ylabel(0(1)4, labsize(medsmall) format(%4.0f) nogrid angle(0)) ///
zlabel(0(0.2)1, labsize(medsmall) format(%4.1f) nogrid) ///
yscale(noline titlegap(1) lwidth(none)) ///
xscale(noline titlegap(1) lwidth(none)) ///
xsize(5) ysize(4) ///
clegend(title("Probability", size(medsmall)) ) ///
text(5.8 0.6 "$ 0.63", color(gs5) placement(e) box size(.4cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
saving(twsa_ext_syrp_1970, replace)



graph combine "twsa_ext_mnp_200" "twsa_ext_syrp_200", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_ext_combined_200, replace) graphregion(color(white))
graph combine "twsa_ext_mnp_985" "twsa_ext_syrp_985", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_ext_combined_985, replace) graphregion(color(white))
graph combine "twsa_ext_mnp_1970" "twsa_ext_syrp_1970", col(2) xsize(8) ysize(4) scheme(s1mono) imargin(0 0 0 0) saving(twsa_ext_combined_1970, replace) graphregion(color(white))
graph export Graph_extended_twsa_1970.tif, replace

graph combine "twsa_ext_combined_200" "twsa_ext_combined_985", col(1) xsize(7.5) ysize(6.25) scheme(s1mono) imargin(0 0 0 0) saving(twsa_ext_combined, replace) plotregion(margin(none))graphregion(color(white))

*graph combine "twsa_combined_200" "twsa_combined_985" "twsa_combined_1970", col(1) xsize(7.5) ysize(9.375) scheme(s1mono) imargin(0 0 0 0) saving(twsa_combined2, replace) plotregion(margin(none))graphregion(color(white))


graph export Graph_extended_twsa.tif, replace
graph export Graph_extended_twsa.eps, replace




cap erase twsa_mnp_200.gph
cap erase twsa_syrp_200.gph
cap erase twsa_mnp_985.gph
cap erase twsa_syrp_985.gph
cap erase twsa_combined_200.gph
cap erase twsa_combined_985.gph
cap erase twsa_mnp_1970.gph
cap erase twsa_syrp_1970.gph
cap erase twsa_combined.gph
cap erase twsa_combined_1970.gph

