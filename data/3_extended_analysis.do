* PART 3:  SECONDARY ANALYSIS USING COMPLETE CASES
*****************************************************************
*****************************************************************
use Analysis_master_data.dta, clear

** Check number of children used for analysis (main paper)
** 662//672//671 at midline and 686//666//684 at endline
drop if hb == .
tabulate treat Rnd

** generate available data and missingness table
bys ScreeningID (treat): egen count = count(ScreeningID)
egen sum_Rnd = sum(Rnd), by(ScreeningID)

generate available = 1 if count == 1 & sum_Rnd == 1
replace  available = 2 if count == 1 & sum_Rnd == 2
replace  available = 3 if count == 1 & sum_Rnd == 3
replace  available = 4 if count == 2 & sum_Rnd == 3
replace  available = 5 if count == 2 & sum_Rnd == 4
replace  available = 6 if count == 2 & sum_Rnd == 5
replace  available = 7 if count == 3 & sum_Rnd == 6
drop sum_Rnd count

#delimit ;
label define AVAILABLE    
     1 "Baseline data only"
     2 "Post-intervention data only"
     3 "Post-follow-up data only"
     4 "Baseline and post-intervention"
     5 "Baseline and post-follow-up"
     6 "Post-intervention and post-follow-up"
     7 "All timepoints (complete cases)"
;
#delimit cr
label values available "AVAILABLE"
tabulate available // complete cases: available==7

** generate anaemia and disability weight variables
drop anaemia
generate anaemia = 1 if hb < 110
recode anaemia . = 0

generate anaemiaCat = 0 if hb >= 110
replace  anaemiaCat = 1 if hb >= 100 & hb < 110
replace  anaemiaCat = 2 if hb >= 70 & hb < 100
replace  anaemiaCat = 3 if hb < 70

generate dw = 0.000 if anaemiaCat == 0
replace  dw = 0.004 if anaemiaCat == 1
replace  dw = 0.052 if anaemiaCat == 2
replace  dw = 0.149 if anaemiaCat == 3

** generate average time-weighted disability weight
bys ScreeningID (Rnd): gen dw4 = ((dw[1]+dw[2])/2)*0.25 if Rnd == 1 & available == 4
bys ScreeningID (Rnd): gen dw5 = ((dw[1]+dw[2])/2)*1.00 if Rnd == 1 & available == 5
bys ScreeningID (Rnd): gen dw6 = ((dw[1]+dw[2])/2)*0.75 if Rnd == 2 & available == 6
bys ScreeningID (Rnd): gen dw7 = ((dw[1]+dw[2])/2)*0.25 + ((dw[2]+dw[3])/2)*0.75 if Rnd == 1 & available == 7

egen dw_all = rowtotal(dw4 dw5 dw6 dw7)
by ScreeningID (Rnd): egen daly = max(dw_all)
replace daly = daly + yld_diarrhoea
drop dw*

** change to negative outcome to ease ICER calculation
replace daly = - daly

save Analysis_complete_data.dta, replace

** Use complete cases: available==7.
use Analysis_complete_data.dta if available == 7 & Rnd == 1, clear

* Total costs from health system perspective
summarize cost_hs if treat == 0
summarize cost_hs if treat == 1
summarize cost_hs if treat == 2

* Total DALYs
summarize daly if treat == 0
summarize daly if treat == 1
summarize daly if treat == 2



** set up bootstrap program

cap program drop boot
cap log close
log using Output_complete_bootstrap_results.log, replace

program define boot

sum cost if treat == 0, meanonly
scalar cost0 = r(mean)
sum cost if treat == 1, meanonly
scalar cost1 = r(mean)
sum cost if treat == 2, meanonly
scalar cost2 = r(mean)

sum cost_hs if treat == 0, meanonly
scalar cost0_hs = r(mean)
sum cost_hs if treat == 1, meanonly
scalar cost1_hs = r(mean)
sum cost_hs if treat == 2, meanonly
scalar cost2_hs = r(mean)

sum daly if treat == 0 & Rnd == 1, meanonly
scalar daly0 = r(mean)
sum daly if treat == 1 & Rnd == 1, meanonly
scalar daly1 = r(mean)
sum daly if treat == 2 & Rnd == 1, meanonly
scalar daly2 = r(mean)

 end


bootstrap cost0=(cost0) cost1=(cost1) cost2=(cost2)    ///
	   cost0_hs=(cost0_hs) cost1_hs=(cost1_hs) cost2_hs=(cost2_hs) ///
          daly0=(daly0) daly1=(daly1) daly2=(daly2)    ///
         icost1=(cost1-cost0) icost2=(cost2-cost0) icost3=(cost2-cost1) ///
      icost1_hs=(cost1_hs-cost0_hs) ///
	  icost2_hs=(cost2_hs-cost0_hs) ///
	  icost3_hs=(cost2_hs-cost1_hs) ///
         idaly1=(daly1-daly0) idaly2=(daly2-daly0) idaly3=(daly2-daly1)  ///
          icer1=((cost1-cost0)/(daly1-daly0))  ///
          icer2=((cost2-cost0)/(daly2-daly0)) ///
          icer3=((cost2-cost1)/(daly2-daly1)) ///
       icer1_hs=((cost1_hs-cost0_hs)/(daly1-daly0))  ///
       icer2_hs=((cost2_hs-cost0_hs)/(daly2-daly0)) ///
       icer3_hs=((cost2_hs-cost1_hs)/(daly2-daly1)), ///
	   reps(2000) strata(treat) seed(12345) ///
	   saving(Analysis_complete_data_bootstrap.dta, replace):boot
estat bootstrap, all		  
cap log close

use Analysis_complete_data_bootstrap.dta, clear


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

egen mcost0_hs = mean(cost0_hs)
egen mcost1_hs = mean(cost1_hs)
egen mcost2_hs = mean(cost2_hs)
egen pcost1_hs = mean(icost1_hs)
egen pcost2_hs = mean(icost2_hs)

save Analysis_complete_data_bootstrap.dta, replace





* PART 4:  SECONDARY ANALYSIS USING COMPLETE CASES AND EXTRAPOLATED
*****************************************************************
*****************************************************************
** Use complete cases: available==7.
use Analysis_complete_data.dta if available == 7, clear

generate dw = 0.000 if anaemiaCat == 0
replace  dw = 0.004 if anaemiaCat == 1
replace  dw = 0.052 if anaemiaCat == 2
replace  dw = 0.149 if anaemiaCat == 3

** generate average time-weighted disability weight
bys ScreeningID (Rnd): gen dw_all = ((dw[1]+dw[2])/2)*0.25 + ((dw[2]+dw[3])/2)*0.75 + ((dw[3]+dw[3])/2)*0.25
drop daly
by ScreeningID (Rnd): egen daly = max(dw_all)
replace daly = daly + yld_diarrhoea
replace daly = - daly
drop dw*

save Analysis_extended_data.dta, replace


* Total costs from health system perspective
summarize cost_hs if treat == 0
summarize cost_hs if treat == 1
summarize cost_hs if treat == 2

* Total DALYs
summarize daly if treat == 0
summarize daly if treat == 1
summarize daly if treat == 2



** set up bootstrap program

cap program drop boot
cap log close
log using Output_extended_bootstrap_results.log, replace

program define boot

sum cost if treat == 0, meanonly
scalar cost0 = r(mean)
sum cost if treat == 1, meanonly
scalar cost1 = r(mean)
sum cost if treat == 2, meanonly
scalar cost2 = r(mean)

sum cost_hs if treat == 0, meanonly
scalar cost0_hs = r(mean)
sum cost_hs if treat == 1, meanonly
scalar cost1_hs = r(mean)
sum cost_hs if treat == 2, meanonly
scalar cost2_hs = r(mean)

sum daly if treat == 0 & Rnd == 1, meanonly
scalar daly0 = r(mean)
sum daly if treat == 1 & Rnd == 1, meanonly
scalar daly1 = r(mean)
sum daly if treat == 2 & Rnd == 1, meanonly
scalar daly2 = r(mean)

 end


bootstrap cost0=(cost0) cost1=(cost1) cost2=(cost2)    ///
	   cost0_hs=(cost0_hs) cost1_hs=(cost1_hs) cost2_hs=(cost2_hs) ///
          daly0=(daly0) daly1=(daly1) daly2=(daly2)    ///
         icost1=(cost1-cost0) icost2=(cost2-cost0) icost3=(cost2-cost1) ///
      icost1_hs=(cost1_hs-cost0_hs) ///
	  icost2_hs=(cost2_hs-cost0_hs) ///
	  icost3_hs=(cost2_hs-cost1_hs) ///
         idaly1=(daly1-daly0) idaly2=(daly2-daly0) idaly3=(daly2-daly1)  ///
          icer1=((cost1-cost0)/(daly1-daly0))  ///
          icer2=((cost2-cost0)/(daly2-daly0)) ///
          icer3=((cost2-cost1)/(daly2-daly1)) ///
       icer1_hs=((cost1_hs-cost0_hs)/(daly1-daly0))  ///
       icer2_hs=((cost2_hs-cost0_hs)/(daly2-daly0)) ///
       icer3_hs=((cost2_hs-cost1_hs)/(daly2-daly1)), ///
	   reps(2000) strata(treat) seed(12345) ///
	   saving(Analysis_extended_data_bootstrap.dta, replace):boot
estat bootstrap, all		  
cap log close

use Analysis_extended_data_bootstrap.dta, clear


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

egen mcost0_hs = mean(cost0_hs)
egen mcost1_hs = mean(cost1_hs)
egen mcost2_hs = mean(cost2_hs)
egen pcost1_hs = mean(icost1_hs)
egen pcost2_hs = mean(icost2_hs)

save Analysis_extended_data_bootstrap.dta, replace



*** Threshold analysis ***************************
** range of values for mnp strategy at 985
postfile contour cost_prgm pce using "thres_ext_mnp_985.dta", replace 
forvalues cost_prgm = 0(0.05)10 {
		use Analysis_complete_data.dta if available == 7 & treat != 2, clear
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 985 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`pce')
}
postclose contour

** range of values for iron strategy at 985
postfile contour cost_prgm pce using "thres_ext_syrp_985.dta", replace 
forvalues cost_prgm = 0(0.05)10 {
		use Analysis_complete_data.dta if available == 7  & treat != 1, clear
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 985 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`pce')
}
postclose contour

** range of values for mnp strategy at 200
postfile contour cost_prgm pce using "thres_ext_mnp_200.dta", replace 
forvalues cost_prgm = 0(0.05)10 {
		use Analysis_complete_data.dta if available == 7 & treat != 2, clear
		replace cost_prgm = `cost_prgm' if treat == 1
		replace cost_intv = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 200 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`pce')
}
postclose contour

** range of values for iron strategy at 200
postfile contour cost_prgm pce using "thres_ext_syrp_200.dta", replace 
forvalues cost_prgm = 0(0.05)10 {
		use Analysis_complete_data.dta if available == 7  & treat != 1, clear
		replace cost_prgm  = `cost_prgm' if treat == 2
		replace cost_intv  = cost_syrp + cost_mnp + cost_prgm
		replace cost_hs = cost_hosp + cost_intv
		gen nmb = daly * 200 - cost_hs
		qui reg nmb treat
		matrix res = r(table)
		local pce = normal(res[1,1]/res[2,1])
		post contour (`cost_prgm') (`pce')
}
postclose contour


use thres_ext_mnp_200.dta, clear
generate comparison = 1
append using thres_ext_syrp_200.dta
replace comparison = 2 if comparison ==.
generate lambda = 200
append using thres_ext_mnp_985.dta
replace comparison = 1 if comparison ==.
append using thres_ext_syrp_985.dta
replace comparison = 2 if comparison ==.
replace lambda = 985 if lambda == .


label define COMPARISON 1 "MNP vs DN" 2 "Iron vs DN"
label values comparison COMPARISON
label variable comparison "Comparison"
label variable pce "Probability cost-effective"
label variable cost_prgm "Program delivery cost"
label variable lambda "Threshold level"
save Analysis_extended_data_thres.dta, replace

cap erase thres_ext_mnp_200.dta
cap erase thres_ext_syrp_200.dta
cap erase thres_ext_mnp_985.dta
cap erase thres_ext_syrp_985.dta



use Analysis_extended_data_thres.dta if lambda == 985, clear

twoway (pci 0 5.65 0.5 5.65, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) (pci 0 3.25 0.5 3.25, lpattern(shortdash) lcolor(gs5) lwidth(vthin))  (pci 0.5 5.65 0.5 0, lpattern(shortdash) lcolor(gs5) lwidth(vthin)) (line pce cost_prgm if comparison == 2, lcolor(dkorange) lwidth(medium) lpattern(solid)) (line pce cost_prgm if comparison == 1, lcolor(emerald) lwidth(medium) lpattern(longdash)), ///
ytitle("Probability cost-effective ({&lambda}=985)", margin(none)) xtitle("Program delivery cost ($)", margin(none)) ///
legend(rows(1) size(small) region(lwidth(none) color(white%0) ) order(4 "Iron supplements" 5 "MNPs") symxsize(6)) ///
yscale(titlegap(2) lwidth(thin)) xscale(titlegap(2) lwidth(thin)) ///
graphregion(lcolor(white) fcolor(white)  margin(medsmall)) ///
plotregion(lwidth(none) fcolor(white) lcolor(gs5)) ///
xsize(16) ysize(13) ///
xlabel(0(1)10, labsize(small) format(%9.0gc) nogrid glcolor(black%3)) ///
ylabel(0(0.1)1,    labsize(small) format(%2.1fc) nogrid glcolor(black%3) angle(0))
graph export Graph_extended_thres.tif, replace
graph export Graph_extended_thres.eps, replace


// xline(5.80,  lpattern(shortdash) lcolor(gs5) lwidth(vthin)) ///
// text(0.98 5.8 "$ 5.80", color(gs5) placement(ne) box size(.7cm)  margin(l+1 r+1 t+1 b+1)  bcolor(white%0) lcolor(white) ) 


