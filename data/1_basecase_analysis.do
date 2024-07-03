* PART 2:  BASE CASE ANALYSIS USING IMPUTED DATA
*****************************************************************
*****************************************************************

use Analysis_master_data.dta, clear

*** log transform ferritin and keep only baseline
gen ferritinL = log(ferritin)
bys ScreeningID (Rnd): replace ferritinL = ferritinL[1] if Rnd != 1
lab var ferritinL "Ferritin level at baseline (log)"
drop ferritin
gen ferritinB = exp(ferritinL)
lab var ferritinB "Ferritin level at baseline"

bys ScreeningID (Rnd): replace anaemia = anaemia[1] if Rnd != 1
rename anaemia anaemia_base
bys ScreeningID (Rnd): replace ironD = ironD[1] if Rnd != 1
rename ironD ironD_base
bys ScreeningID (Rnd): replace ironDA = ironDA[1] if Rnd != 1
rename ironDA ironDA_base

lab var anaemia_base "Anaemia at baseline"
lab var ironD_base   "Iron deficiency at baseline"
lab var ironDA_base  "Iron deficiency anaemia at baseline"

*** check number of missing values
misstable summarize

*** mean imputation for baseline variables (just a few % missing)
summ Qa9_cat
replace Qa9_cat = round(r(mean)) if Qa9_cat ==.

summ wealth
replace wealth = r(mean) if wealth ==.

bys ScreeningID (fci_total): egen fci_total_mean = mean(fci_total) 
replace fci_total = fci_total_mean if fci_total == . & fci_total_mean != .
drop fci_total_mean

summ fci_total
replace fci_total = r(mean) if fci_total == .


*** Reshape data to wide for "just another variable" imputation
reshape wide hb, i(ScreeningID) j(Rnd)
egen hb_count = rownonmiss(hb1 hb2 hb3)
gen complete = 1 if hb_count == 3
replace complete = 0 if complete == .

lab var fci_total "Family care indicator"
lab var compliance "70% Compliance"
lab var anaemia_base "Baseline anaemia"
lab var ironD_base "Baseline iron deficiency"
lab var ironDA_base "Baseline iron deficiency anaemia"
lab define COMPLETE 0 "Missing" 1 "Complete"
lab values complete COMPLETE
tab complete

cap log close
log using Output_basechar_missing.log, replace

*** Table 3. Baseline Characteristics of Participants by treatment group and data availability
table1_mc, by(treat) vars(female bin \ wealth_quintile cat \ fci_total contn %4.1f \ Qa9_cat cat \ union cat \ anaemia_base bin \ ironD_base bin \ ironDA_base bin \ compliance bin)  onecol slashN
 
table1_mc if treat==0, by(complete) vars(female bin \ wealth_quintile cat \ fci_total contn %4.1f \ Qa9_cat cat \ union cat \ anaemia_base bin \ ironD_base bin \ ironDA_base bin \ compliance bin)  onecol slashN
 
table1_mc if treat==1, by(complete) vars(female bin \ wealth_quintile cat \ fci_total contn %4.1f \ Qa9_cat cat \ union cat \ anaemia_base bin \ ironD_base bin \ ironDA_base bin \ compliance bin)  onecol slashN

table1_mc if treat==2, by(complete) vars(female bin \ wealth_quintile cat \ fci_total contn %4.1f \ Qa9_cat cat \ union cat \ anaemia_base bin \ ironD_base bin \ ironDA_base bin \ compliance bin)  onecol slashN

cap log close

drop if hb1 == .
*drop if hb_count < 2
drop hb_count complete
table1_mc, by(treat) vars(hb1 contn \ hb2 contn \ hb3 contn)
di 422/1063   // = 39.8%. Use about 40 imputations

*** Set up for multiple imputation
mi set flong
mi register imputed hb* ferritinL*         // requires imputation
mi register regular union female fci_total Qa9_cat compliance // not require imputation

*** Implement multiple imputation
mi impute chained (regress) hb2 hb3 ferritinL = i.union i.female fci i.Qa9_cat i.compliance, add(40) by(treat) rseed(12345)
// MI by chained equations using linear regression
// Imputing hb at each Rnd, by treatment group
// Using variables union, female, fci and mEduc
// Perform 40 imputations, stratified by treatment group.
// For reproducibility, set seed 12345

mi passive: gen ferritin = exp(ferritinL)   
// reverse log transformation of ferritin

*** Reshape data back to long format for pooling and analysis
mi reshape long hb, i(ScreeningID) j(Rnd)

*** Generate variables needed in imputed data
mi passive: gen anaemia = 1 if hb < 110
mi passive: replace anaemia = 0 if hb >= 110
label variable anaemia "Anaemia status (binary)"
mi passive: replace anaemia = -anaemia

mi passive: generate anaemiaCat = 0 if hb >= 110
mi passive: replace anaemiaCat = 1 if hb >= 100 & hb < 110
mi passive: replace anaemiaCat = 2 if hb >= 70 & hb < 100
mi passive: replace anaemiaCat = 3 if hb < 70
label variable anaemiaCat "Anaemia severity (categorical)"

*** Generate time weighted disability weights
sort ScreeningID Rnd

* base case
mi passive: generate dw = 0.000 if anaemiaCat == 0
mi passive: replace dw = 0.004 if anaemiaCat == 1
mi passive: replace dw = 0.052 if anaemiaCat == 2
mi passive: replace dw = 0.149 if anaemiaCat == 3
mi passive: by ScreeningID (Rnd): gen dw_est = ((dw[1]+dw[2])/2)*0.25 + ((dw[2]+dw[3])/2)*0.75 if Rnd == 1  // 3 mnth = 0.25, 9 mnth = 0.75
mi passive: by ScreeningID (Rnd): egen daly = max(dw_est)
mi passive: replace daly = daly + yld_diarrhoea
mi passive: replace daly = -daly

* low values
mi passive: generate dw_low = 0.000 if anaemiaCat == 0
mi passive: replace  dw_low = 0.001 if anaemiaCat == 1
mi passive: replace  dw_low = 0.034 if anaemiaCat == 2
mi passive: replace  dw_low = 0.101 if anaemiaCat == 3
mi passive: by ScreeningID (Rnd): gen dw_est_low = ((dw_low[1]+dw_low[2])/2)*0.25 + ((dw_low[2]+dw_low[3])/2)*0.75 if Rnd == 1
mi passive: by ScreeningID (Rnd): egen daly_low = max(dw_est_low)
mi passive: replace daly_low = daly_low + yld_diarrhoea
mi passive: replace daly_low = - daly_low

* high values
mi passive: generate dw_high = 0.000 if anaemiaCat == 0
mi passive: replace  dw_high = 0.008 if anaemiaCat == 1
mi passive: replace  dw_high = 0.076 if anaemiaCat == 2
mi passive: replace  dw_high = 0.209 if anaemiaCat == 3
mi passive: by ScreeningID (Rnd): gen dw_est_high = ((dw_high[1]+dw_high[2])/2)*0.25 + ((dw_high[2]+dw_high[3])/2)*0.75 if Rnd == 1
mi passive: by ScreeningID (Rnd): egen daly_high = max(dw_est_high)
mi passive: replace daly_high = daly_high + yld_diarrhoea
mi passive: replace daly_high = - daly_high

format daly* %9.3f
drop  dw* anaemiaCat ferritinL

save Analysis_imputed_data.dta, replace



***************************************************************************
*** Output: Resource use by treatment arm
cap log close
log using Output_resource_use.log, replace

use Analysis_imputed_data.dta, clear
mi extract 0
keep if Rnd == 1
collapse (sum) num_OPD_visit num_IPD_bed num_IPD_diarr, by(treat)

use Analysis_imputed_data.dta, clear
mi extract 0
sort ScreeningID Rnd
tabstat num_OPD_visit if Rnd == 1, by(treat) stat(mean sd n)
tabstat num_IPD_bed if Rnd == 1, by(treat) stat(mean sd n)
tabstat num_IPD_diarr if Rnd == 1, by(treat) stat(mean sd n)

cap log close


*** Output: Costs by treatment arm
cap log close
log using Output_costs_perspectives.log, replace

use Analysis_imputed_data.dta, clear
mi extract 0
bys treat: summ cost_mnp  if Rnd == 1
bys treat: summ cost_syrp if Rnd == 1
bys treat: summ cost_prgm if Rnd == 1

* Cost of OPD visit due to diarrhoea
ttest cost_OPD_visit if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost_OPD_visit if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost_OPD_visit if treat !=0 & Rnd == 1, by(treat) reverse

* Cost of OPD treatment of diarrhoea
ttest cost_OPD_tx if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost_OPD_tx if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost_OPD_tx if treat !=0 & Rnd == 1, by(treat) reverse

* Cost of IPD be days due to diarrhoea
ttest cost_IPD_bed if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost_IPD_bed if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost_IPD_bed if treat !=0 & Rnd == 1, by(treat) reverse

* Cost of IPD treatment of diarrhoea
ttest cost_IPD_tx if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost_IPD_tx if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost_IPD_tx if treat !=0 & Rnd == 1, by(treat) reverse

* Total costs from health system perspective
ttest cost_hs if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost_hs if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost_hs if treat !=0 & Rnd == 1, by(treat) reverse

* Costs out-of-pocket for diarrhoea treatment
ttest cost_oop if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost_oop if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost_oop if treat !=0 & Rnd == 1, by(treat) reverse

* Costs due to lost productivity by caregivers
ttest cost_prod if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost_prod if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost_prod if treat !=0 & Rnd == 1, by(treat) reverse

* Total costs from societal perspective
ttest cost if treat !=2 & Rnd == 1, by(treat) reverse
ttest cost if treat !=1 & Rnd == 1, by(treat) reverse
ttest cost if treat !=0 & Rnd == 1, by(treat) reverse

cap log close

***************************************************************************

use Analysis_imputed_data.dta, clear
set cformat %9.4f


summarize daly if treat == 0 & Rnd == 1
summarize daly if treat == 1 & Rnd == 1
summarize daly if treat == 2 & Rnd == 1



*** Bootstrap + Imputation: Health system perspective

* data in ICE format (requires the stata package ice)


mi export ice
keep if Rnd == 1
keep ScreeningID treat cost cost_hs _mi _mj daly
save Analysis_imputed_data_ice.dta, replace



mim: sum daly if treat == 0
matrix define d0 = e(MIM_Q)
scalar daly0 = d0[1,1]

mim: mean daly if treat == 1
matrix define d1 = e(MIM_Q)
scalar daly1 = d1[1,1]

mim: mean daly if treat == 2
matrix define d2 = e(MIM_Q)
scalar daly2 = d2[1,1]

cap prog drop boot
program define boot

mim: mean cost_hs if treat == 0
matrix define c0_HS = e(MIM_Q)
scalar cost0_HS = c0_HS[1,1]

mim: mean cost_hs if treat == 1
matrix define c1_HS = e(MIM_Q)
scalar cost1_HS = c1_HS[1,1]

mim: mean cost_hs if treat == 2
matrix define c2_HS = e(MIM_Q)
scalar cost2_HS = c2_HS[1,1]

mim: mean daly if treat == 0
matrix define d0 = e(MIM_Q)
scalar daly0 = d0[1,1]

mim: mean daly if treat == 1
matrix define d1 = e(MIM_Q)
scalar daly1 = d1[1,1]

mim: mean daly if treat == 2
matrix define d2 = e(MIM_Q)
scalar daly2 = d2[1,1]

end

cap log close
log using Output_imputed_bootstrap_results_hs.log, replace
bootstrap cost0_HS=(cost0_HS) cost1_HS=(cost1_HS) cost2_HS=(cost2_HS) ///
            daly0=(daly0) daly1=(daly1) daly2=(daly2) ///
        icost1_HS=(cost1_HS-cost0_HS) ///
	    icost2_HS=(cost2_HS-cost0_HS) ///
	    icost3_HS=(cost2_HS-cost1_HS) ///
           idaly1=(daly1-daly0) idaly2=(daly2-daly0) idaly3=(daly2-daly1) ///
         icer1_HS=((cost1_HS-cost0_HS)/(daly1-daly0))  ///
         icer2_HS=((cost2_HS-cost0_HS)/(daly2-daly0)) ///
         icer3_HS=((cost2_HS-cost1_HS)/(daly2-daly1)), ///
	   reps(2000) strata(treat) seed(12345) cluster(_mi) idcluster(newvar)  ///
	   saving(Analysis_imputed_data_bootstrap_hs.dta, replace):boot
estat bootstrap, all		  
cap log close

use Analysis_imputed_data_bootstrap_hs.dta, clear

** generate difference in costs
egen mcost0_HS = mean(cost0_HS)
egen mcost1_HS = mean(cost1_HS)
egen mcost2_HS = mean(cost2_HS)
egen mdaly0 = mean(daly0)
egen mdaly1 = mean(daly1)
egen mdaly2 = mean(daly2)
egen pcost1_HS = mean(icost1_HS)
egen pcost2_HS = mean(icost2_HS)
egen pdaly1 = mean(idaly1)
egen pdaly2 = mean(idaly2)
gen  icost0 = cost0 - cost0
gen  idaly0 = daly0 - daly0

save Analysis_imputed_data_bootstrap_hs.dta, replace



*** Bootstrap + Imputation: Societal perspective
* data in ICE format (requires the stata package ice)
use Analysis_imputed_data.dta, clear
mi export ice
keep if Rnd == 1
keep ScreeningID treat cost cost_hs _mi _mj daly
save Analysis_imputed_data_ice.dta, replace

cap prog drop boot
program define boot

mim: mean cost if treat == 0
matrix define c0 = e(MIM_Q)
scalar cost0 = c0[1,1]

mim: mean cost if treat == 1
matrix define c1 = e(MIM_Q)
scalar cost1 = c1[1,1]

mim: mean cost if treat == 2
matrix define c2 = e(MIM_Q)
scalar cost2 = c2[1,1]

mim: mean daly if treat == 0
matrix define d0 = e(MIM_Q)
scalar daly0 = d0[1,1]

mim: mean daly if treat == 1
matrix define d1 = e(MIM_Q)
scalar daly1 = d1[1,1]

mim: mean daly if treat == 2
matrix define d2 = e(MIM_Q)
scalar daly2 = d2[1,1]

end

cap log close
log using Output_imputed_bootstrap_results.log, replace
bootstrap cost0=(cost0) cost1=(cost1) cost2=(cost2) ///
          daly0=(daly0) daly1=(daly1) daly2=(daly2) ///
         icost1=(cost1-cost0) icost2=(cost2-cost0) icost3=(cost2-cost1) ///
         idaly1=(daly1-daly0) idaly2=(daly2-daly0) idaly3=(daly2-daly1) ///
          icer1=((cost1-cost0)/(daly1-daly0))  ///
          icer2=((cost2-cost0)/(daly2-daly0)) ///
          icer3=((cost2-cost1)/(daly2-daly1)), ///
	   reps(2000) strata(treat) seed(12345) cluster(_mi) ///
	   saving(Analysis_imputed_data_bootstrap.dta, replace):boot
estat bootstrap, all		  
cap log close

use Analysis_imputed_data_bootstrap.dta, clear

** generate difference in costs
egen mcost0 = mean(cost0)
egen mcost1 = mean(cost1)
egen mcost2 = mean(cost2)
egen mdaly0 = mean(daly0)
egen mdaly1 = mean(daly1)
egen mdaly2 = mean(daly2)
egen pcost1 = mean(icost1)
egen pcost2 = mean(icost2)
egen pdaly1 = mean(idaly1)
egen pdaly2 = mean(idaly2)
gen  icost0 = cost0 - cost0
gen  idaly0 = daly0 - daly0

save Analysis_imputed_data_bootstrap.dta, replace



*** COST-EFFECTIVENESS PLANE ***************************
use Analysis_imputed_data_bootstrap_hs.dta, clear
set scheme s1color

range x200 -0.001 0.005 2000
gen y200 = 200*x200
range x985 -0.001 0.005 2000
gen y985 = 985*x985
range x1970 -0.001 0.0045 2000
gen y1970 = 1970*x1970

* health system perspective for syrp and mnp vs dn

twoway (scatter icost1_HS idaly1, mcolor(emerald%90) msize(medium) msymbol(circle) mfcolor(white%80)) ///
(scatter icost2_HS idaly2, mcolor(dkorange) msize(medium) msymbol(circle) mfcolor(white%80)) ///
(scatter pcost1_HS pdaly1, mcolor(emerald%90) msize(medium) msymbol(circle)) ///
(scatter pcost2_HS pdaly2, mcolor(dkorange) msize(medium) msymbol(circle)) ///
(line y200  x200,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
(line y985  x985,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
(line y1970 x1970, lpattern(shortdash) lcolor(gs5) lwidth(vthin)), ///
ytitle("Incremental costs ($)", margin(none)) xtitle("DALYs averted", margin(none)) ///
plotregion(lwidth(none) fcolor(white)) ///
legend(rows(1) size(small)  region(lwidth(none) color(white%0)) order(3 "MNPs vs No intervention" 4 "Iron supplements vs No intervention" ) position(5) ) ///
graphregion(lcolor(white) fcolor(white)  margin(small)) ///
yline(0, lc(black) lwidth(medthin)) xline(0, lc(black) lwidth(medthin)) ///
xsize(16) ysize(13) ///
ylabel(-2(2)8,             labsize(small) ticks nogrid glcolor(black%3) format(%9.0g) angle(0)) ///
xlabel(-0.001(0.001)0.006, labsize(small) ticks nogrid glcolor(black%3) ) ///
yscale(noline titlegap(2) lwidth(none)) xscale(noline titlegap(2) lwidth(none)) ///
text(1 0.005 "{&lambda} = $ 200", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(5 0.005 "{&lambda} = $ 985", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(9 0.0045 "{&lambda} = $ 1970", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) 
graph export Graph_imputed_plane_both_vs_dn_hs.tif, replace
graph export Graph_imputed_plane_both_vs_dn_hs.eps, replace



* health system perspective for syrp and mnp vs dn
twoway (scatter icost3_HS idaly3, mcolor(gs5%90) msize(medium) msymbol(circle) mfcolor(white%80)), ///
ytitle("Incremental costs ($)", margin(none)) xtitle("DALYs averted", margin(none)) ///
plotregion(lwidth(vvthin) fcolor(white)) ///
legend(rows(1) size(small)  region(lwidth(none) color(white%0)) order(1 "Iron supplements vs MNPs") position(5) ) ///
graphregion(lcolor(white) fcolor(white)  margin(small)) ///
yline(0, lc(black) lwidth(medthin)) xline(0, lc(black) lwidth(medthin)) ///
xsize(16) ysize(13) ///
ylabel(-2(2)8, labsize(small) ticks nogrid  glcolor(black%3) format(%9.0g) angle(0)) ///
xlabel(-0.001(0.001)0.006,         labsize(small) ticks nogrid  glcolor(black%3) ) ///
yscale(noline titlegap(2) lwidth(none)) xscale(noline titlegap(2) lwidth(none)) 
graph export Graph_imputed_plane_syrup_vs_mnp_hs.tif, replace
graph export Graph_imputed_plane_syrup_vs_mnp_hs.eps, replace






use Analysis_imputed_data_bootstrap.dta, clear
set scheme s1color

range x200 -0.001 0.005 2000
gen y200 = 200*x200
range x985 -0.001 0.005 2000
gen y985 = 985*x985
range x1970 -0.001 0.0045 2000
gen y1970 = 1970*x1970


* societal perspective for syrp and mnp vs dn
twoway (scatter icost1 idaly1, mcolor(emerald%90) msize(medium) msymbol(circle) mfcolor(white%80)) ///
(scatter icost2 idaly2, mcolor(dkorange) msize(medium) msymbol(circle) mfcolor(white%80)) ///
(scatter pcost1 pdaly1, mcolor(emerald%90) msize(medium) msymbol(circle)) ///
(scatter pcost2 pdaly2, mcolor(dkorange) msize(medium) msymbol(circle)) ///
(line y200  x200,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
(line y985  x985,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
(line y1970 x1970, lpattern(shortdash) lcolor(gs5) lwidth(vthin)), ///
ytitle("Incremental costs ($)", margin(none)) xtitle("DALYs averted", margin(none)) ///
plotregion(lwidth(vvthin) fcolor(white)) ///
legend(rows(1) size(small)  region(lwidth(none) color(white%0)) order(3 "MNPs vs No intervention" 4 "Iron supplements vs No intervention" ) position(5) ) ///
graphregion(lcolor(white) fcolor(white)  margin(small)) ///
yline(0, lc(black) lwidth(medthin)) xline(0, lc(black) lwidth(medthin)) ///
xsize(16) ysize(13) ///
ylabel(-2(2)8, labsize(small) ticks nogrid  glcolor(black%3) format(%9.0g) angle(0)) ///
xlabel(-0.001(0.001)0.006,         labsize(small) ticks nogrid  glcolor(black%3) ) ///
yscale(noline titlegap(2) lwidth(none)) xscale(noline titlegap(2) lwidth(none)) ///
text(1 0.005 "{&lambda} = $ 200", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(5 0.005 "{&lambda} = $ 985", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(9 0.0045 "{&lambda} = GDP per capita", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) 
graph export Graph_imputed_plane_both_vs_dn.tif, replace
graph export Graph_imputed_plane_both_vs_dn.eps, replace


* societal perspective for syrp and mnp vs dn
twoway (scatter icost3 idaly3, mcolor(gs5%90) msize(medium) msymbol(circle) mfcolor(white%80)), ///
ytitle("Incremental costs ($)", margin(none)) xtitle("DALYs averted", margin(none)) ///
plotregion(lwidth(vvthin) fcolor(white)) ///
legend(rows(1) size(small)  region(lwidth(none) color(white%0)) order(1 "Iron supplements vs MNPs") position(5) ) ///
graphregion(lcolor(white) fcolor(white)  margin(small)) ///
yline(0, lc(black) lwidth(medthin)) xline(0, lc(black) lwidth(medthin)) ///
xsize(16) ysize(13) ///
ylabel(-2(2)8, labsize(small) ticks grid  glcolor(black%3) format(%9.0g) angle(0)) ///
xlabel(-0.001(0.001)0.006,         labsize(small) ticks nogrid  glcolor(black%3) ) ///
yscale(noline titlegap(2) lwidth(none)) xscale(noline titlegap(2) lwidth(none)) 
graph export Graph_imputed_plane_syrup_vs_mnp.tif, replace
graph export Graph_imputed_plane_syrup_vs_mnp.eps, replace


*** COST-EFFECTIVENESS ACCEPTABILITY CURVES ***************************
use Analysis_imputed_data_bootstrap_hs.dta, clear

postfile ceac wtp prob0_HS prob1_HS prob2_HS using "Analysis_imputed_data_ceac_hs.dta", replace
forvalues wtp = 0(1)3000 {
	count if ((daly0*`wtp'- cost0_HS) > (daly1*`wtp'- cost1_HS)) &  ((daly0*`wtp'- cost0_HS) > (daly2*`wtp'- cost2_HS)) 
	local p0_HS = `r(N)' / _N 
	count if ((daly1*`wtp'- cost1_HS) > (daly0*`wtp'- cost0_HS)) &  ((daly1*`wtp'- cost1_HS) > (daly2*`wtp'- cost2_HS)) 
	local p1_HS = `r(N)' / _N 
	count if ((daly2*`wtp'- cost2_HS) > (daly0*`wtp'- cost0_HS)) &  ((daly2*`wtp'- cost2_HS) > (daly1*`wtp'- cost1_HS)) 
	local p2_HS = `r(N)' / _N 
	post ceac (`wtp') (`p0_HS') (`p1_HS') (`p2_HS')
	}
postclose ceac //Closing postfile

use Analysis_imputed_data_ceac_hs.dta, clear
label var wtp   "Willingness to pay threshold"
label var prob0_HS "Probability doing nothing is preferred (Health system)"
label var prob1_HS "Probability fortification is preferred (Health system)"
label var prob2_HS "Probability supplementation is preferred (Health system)"
save Analysis_imputed_data_ceac_hs.dta, replace


use Analysis_imputed_data_ceac_hs.dta, clear

** Health system perspective for iron vs placebo & mnp vs placebo
* Coloured
twoway (line prob0_HS wtp, lcolor(gs0) lwidth(medium) lpattern(solid))  (line prob2_HS wtp, lcolor(dkorange) lwidth(medthick) lpattern(shortdash)) (line prob1_HS wtp, lcolor(emerald) lwidth(medthick) lpattern(longdash)), ///
xline(200,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
xline(985,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
xline(1970, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
ytitle("Probability cost-effective", margin(none)) xtitle("Cost per DALY averted threshold ($)", margin(none)) ///
legend(rows(1) pos(5) size(small) region(lwidth(none) color(white%0) ) order(1 "No intervention" 3 "MNPs" 2 "Iron supplements") symxsize(9)) ///
yscale(titlegap(2) lwidth(thin)) xscale(titlegap(2) lwidth(thin)) ///
graphregion(lcolor(white) fcolor(white)  margin(medsmall)) ///
plotregion(lwidth(none) fcolor(white) lcolor(gs5)) ///
xsize(48) ysize(39) ///
xlabel(0(500)3000, labsize(3) format(%9.0gc) nogrid glcolor(black%3)) ///
ylabel(0(0.2)1,    labsize(3) format(%2.1fc) nogrid glcolor(black%3) angle(0)) ///
text(0.49 200 "$ 200", color(gs5) placement(ne) box size(3cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(0.49 985 "$ 985", color(gs5) placement(ne) box size(3cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(0.49 1970 "$ 1970", color(gs5) placement(ne) box size(3cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) )
graph export figure1_color.png, replace
graph export figure1_color.eps, replace

* Black and white
twoway (line prob0_HS wtp, lcolor(black) lwidth(medium) lpattern(solid))  (line prob2_HS wtp, lcolor(black) lwidth(medthick) lpattern(shortdash)) (line prob1_HS wtp, lcolor(black) lwidth(medthick) lpattern(longdash)), ///
xline(200,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
xline(985,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
xline(1970, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
ytitle("Probability cost-effective", margin(none)) xtitle("Cost per DALY averted threshold ($)", margin(none)) ///
legend(rows(1) pos(5) size(small) region(lwidth(none) color(white%0) ) order(1 "No intervention" 3 "MNPs" 2 "Iron supplements") symxsize(9)) ///
yscale(titlegap(2) lwidth(thin)) xscale(titlegap(2) lwidth(thin)) ///
graphregion(lcolor(white) fcolor(white)  margin(medsmall)) ///
plotregion(lwidth(none) fcolor(white) lcolor(gs5)) ///
xsize(48) ysize(39) ///
xlabel(0(500)3000, labsize(3) format(%9.0gc) nogrid glcolor(black%3)) ///
ylabel(0(0.2)1,    labsize(3) format(%2.1fc) nogrid glcolor(black%3) angle(0)) ///
text(0.49 200 "$ 200", color(gs5) placement(ne) box size(3cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(0.49 985 "$ 985", color(gs5) placement(ne) box size(3cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(0.49 1970 "$ 1970", color(gs5) placement(ne) box size(3cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) )
graph export figure1_blackwhite.png, replace
graph export figure1_blackwhite.eps, replace


use Analysis_imputed_data_bootstrap.dta, clear

postfile ceac wtp prob0 prob1 prob2 using "Analysis_imputed_data_ceac.dta", replace
forvalues wtp = 0(1)3000 {
	count if ((daly0*`wtp'- cost0) > (daly1*`wtp'- cost1)) &  ((daly0*`wtp'- cost0) > (daly2*`wtp'- cost2)) 
	local p0 = `r(N)' / _N 
	count if ((daly1*`wtp'- cost1) > (daly0*`wtp'- cost0)) &  ((daly1*`wtp'- cost1) > (daly2*`wtp'- cost2)) 
	local p1 = `r(N)' / _N 
	count if ((daly2*`wtp'- cost2) > (daly0*`wtp'- cost0)) &  ((daly2*`wtp'- cost2) > (daly1*`wtp'- cost1)) 
	local p2 = `r(N)' / _N 
	post ceac (`wtp') (`p0') (`p1') (`p2')
	}
postclose ceac //Closing postfile

use Analysis_imputed_data_ceac.dta, clear
label var wtp   "Willingness to pay threshold"
label var prob0 "Probability doing nothing is preferred"
label var prob1 "Probability fortification is preferred"
label var prob2 "Probability supplementation is preferred"
save Analysis_imputed_data_ceac.dta, replace

use Analysis_imputed_data_ceac.dta, clear


** societal perspective for iron vs placebo & mnp vs placebo
twoway (line prob0 wtp, lcolor(gs0) lwidth(medium) lpattern(solid))  (line prob2 wtp, lcolor(dkorange) lwidth(medium) lpattern(solid)) (line prob1 wtp, lcolor(emerald) lwidth(medium) lpattern(longdash)), ///
xline(200,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
xline(985,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
xline(1970, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
ytitle("Probability cost-effective", margin(none)) xtitle("WTP per DALY averted ($)", margin(none)) ///
legend(rows(1) pos(6) size(small) region(lwidth(none) color(white%0) ) order(1 "No intervention" 3 "MNPs" 2 "Iron supplements") symxsize(6)) ///
yscale(titlegap(2) lwidth(thin)) xscale(titlegap(2) lwidth(thin)) ///
graphregion(lcolor(white) fcolor(white)  margin(medsmall)) ///
plotregion(lwidth(none) fcolor(white) lcolor(gs5)) ///
xsize(16) ysize(13) ///
xlabel(0(500)3000, labsize(small) format(%9.0gc) nogrid glcolor(black%3)) ///
ylabel(0(0.2)1,    labsize(small) format(%2.1fc) nogrid glcolor(black%3)) ///
text(0.8 200 "{&lambda} = $ 200", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(0.8 985 "{&lambda} = $ 985", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) ///
text(0.8 1970 "{&lambda} = GDP per capita", color(gs5) placement(e) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) )
graph export Graph_imputed_ceac.tif, replace
graph export Graph_imputed_ceac.eps, replace




