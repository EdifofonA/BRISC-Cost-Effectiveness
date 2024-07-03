clear all
clear matrix


* PART 1:  CREATE MASTER DATASET (brisc.dta)
*****************************************************************
****************************************************************



***** A. INPUTS / UNIT COSTS / INFLATION / CURRENCY*****
****************************************************************

*** Conversion rates
scalar xch10 =  69.6  // USD to Taka in 2010
scalar xch15 =  77.9  // USD to Taka in 2015
scalar xch17 =  80.4  // USD to Taka in 2017
scalar xch20 =  84.9  // USD to Taka in 2020
scalar def10 = 131.4  // USD to Taka in 2010
scalar def13 = 164.3  // USD to Taka in 2010

scalar def15 = 183.8  // GDP deflator for Bangladesh  in 2015
scalar def17 = 208.4  // GDP deflator for Bangladesh  in 2015
scalar def20 = 242.6  // GDP deflator for Bangladesh  in 2015
scalar ppp10 =  22.2  // International Dollar to Taka in 2008

*** health service costs from WHO-CHOICE 2010 I$ to 2020 US$
scalar unit_bed        = 14.68 * ppp10 * (def20/def10) / xch20 // Base 
scalar unit_visit      =  3.50 * ppp10 * (def20/def10) / xch20 // Base 
scalar unit_bed_low    = unit_bed   / 14.68 *  5.43            // Low 
scalar unit_visit_low  = unit_visit /  3.50 *  0.74            // Low 
scalar unit_bed_high   = unit_bed   / 14.68 * 30.42            // High 
scalar unit_visit_high = unit_visit /  3.50 * 10.32            // High 
   
** cost of medications for treating diarrhoea
scalar ors     =   4.50/xch20  //Cost of ORS sachet (BDT)
scalar ringers =  71.30/xch20  //Cost of Ringers lactate (BDT)
scalar zinc    =   1.75/xch20  //Cost of Zinc tablets (BDT)

** cost of intervention and programme management
scalar unit_prgm  =  4.50 * 74.15 * (def20/def13) / xch20
scalar unit_syrp  =         40.25 * (def20/def15) / xch20
scalar unit_mnp   =          1.30 * (def20/def17) / xch20

scalar dw_diarr_severe     = 0.247  // 0.247 (0.164-0.348) from GBD 2019 
scalar dw_diarr_moderate   = 0.188  // 0.188 (0.125-0.264) from GBD 2019

scalar dur_diarr_severe    = 8.4  // 0.247 (0.164-0.348) from GBD 2019 
scalar dur_diarr_moderate  = 6.4  // 0.188 (0.125-0.264) from GBD 2019


// Lamberti LM,  Fischer  Walker CL,  Black  RE.  Systematic  review  of  diarrhea duration and severity  in children  and adults  in low-  and middle-income  countries.  BMC  Public  Health 2012;  12:  276


***** B. CLEAN INPATIENT COSTS FROM SAE_HE DATASET
*****************************************************************

*** Create file of children and treatment groups (to delete later)
use BaseChar_HE.dta, clear
keep ScreeningID treat
save randomization.dta, replace

*** Use hospital stays where Acute Diarrhoea (AWD) was diagnosis
use SAE_HE.dta, clear
gen num_diarrhoea      = 1 if strmatch(Diagnosis , "*AWD*")
keep if num_diarrhoea == 1

*** Number of nights is one less than number of days
replace duration = duration - 1

*** correct errors in duration of illness in records 63 & 88
replace duration = 1 in 63
replace duration = 1 in 88

*** Use mean imputation for 6 missing duration of stay values
su duration
replace duration   = round(r(mean)) if duration == .

*** Generate total number of bed days per child
egen num_IPD_bed   = total(duration),      by(ScreeningID)
egen num_IPD_diarr = total(num_diarrhoea), by(ScreeningID)
keep ScreeningID num_IPD_bed num_IPD_diarr
duplicates drop

*** Generate cost of inpatient stay = bed days + severe diarrhoea tx
gen cost_IPD_bed = num_IPD_bed   * unit_bed
gen cost_IPD_tx  = num_IPD_diarr * (ringers*2 + ors*3 + zinc*10) 
gen cost_IPD     = cost_IPD_bed + cost_IPD_tx

* Sensitivity analysis: Vary treatment cost by 50% in either direction
gen cost_IPD_tx_low  = cost_IPD_bed + num_IPD_diarr * (ringers*2 + ors*3 + zinc*10)/1.5 
gen cost_IPD_tx_high = cost_IPD_bed + num_IPD_diarr * (ringers*2 + ors*3 + zinc*10)*1.5
* Sensitivity analysis: Vary bed day cost using 95% uncertainty intervals
gen cost_IPD_bed_low  = (num_IPD_bed*unit_bed_low)  + cost_IPD_tx
gen cost_IPD_bed_high = (num_IPD_bed*unit_bed_high)  + cost_IPD_tx

*** Merge with randomization file to get complete 3300 observations
merge 1:1 ScreeningID using randomization.dta
foreach x of varlist num* cost*{
  replace `x' = 0 if missing(`x') //replace missing values with 0
}
drop _merge
sort ScreeningID
order ScreeningID treat

gen yld_diarrhoea = (num_IPD_diarr*dur_diarr_severe/365.25)*dw_diarr_severe
save inpatient.dta, replace



***** C. CLEAN OUTPATIENT COSTS FROM WEEKLYMORB_HE/MONTHLYMORB_HE DATA
*****************************************************************

*** FOR INTERVENTION PHASE ************
*** Use data with with Diarrhoea as reason for last hospital visit
use WeeklyMorb_HE.dta, clear
keep if Diarrhoea == 1

*** Generate total number of hospital visits per child
egen num_OPD_visit = total(Diarrhoea), by(ScreeningID)
gen  num_OPD_diarr = num_OPD_visit

*** Generate cost of outpatient visit = visits + diarrhoea tx
gen cost_OPD_visit = num_OPD_visit * unit_visit
gen cost_OPD_tx    = num_OPD_diarr * (ors*2 + zinc*10)
gen cost_OPD       = cost_OPD_visit + cost_OPD_tx

* Sensitivity analysis: Vary treatment cost by 50% in either direction
gen cost_OPD_tx_low  = num_OPD_diarr * (ors*2 + zinc*10)/1.5
gen cost_OPD_tx_high = num_OPD_diarr * (ors*2 + zinc*10)*1.5

* Sensitivity analysis: Vary visit cost using 95% uncertainty intervals
gen cost_OPD_visit_low  = (num_OPD_visit * unit_visit_low)  + cost_OPD_tx
gen cost_OPD_visit_high = (num_OPD_visit * unit_visit_high) + cost_OPD_tx

keep ScreeningID num* cost*
duplicates drop
save weekly.dta, replace


*** FOR FOLLOW-UP PHASE ************
*** Use data with with Diarrhoea as reason for last hospital visit
use MonthlyMorb_HE.dta, clear
keep if Diarrhoea == 1

*** Generate total number of hospital visits
egen num_OPD_visit = total(Diarrhoea), by(ScreeningID)
replace num_OPD_visit = num_OPD_visit * 2.167 
// rescaled from 14 days to 1 month because data collected at this stage was based on 14 days recall
gen  num_OPD_diarr = num_OPD_visit

gen cost_OPD_visit      = num_OPD_visit * unit_visit
gen cost_OPD_tx      = num_OPD_diarr * (ors*2 + zinc*10)
gen cost_OPD         = cost_OPD_visit + cost_OPD_tx

* Sensitivity analysis: Vary treatment cost by 50% in either direction
gen cost_OPD_tx_low  = num_OPD_diarr * (ors*2 + zinc*10)/1.5
gen cost_OPD_tx_high = num_OPD_diarr * (ors*2 + zinc*10)*1.5

* Sensitivity analysis: Vary visit cost using 95% uncertainty intervals
gen cost_OPD_visit_low  = (num_OPD_visit * unit_visit_low)  + cost_OPD_tx
gen cost_OPD_visit_high = (num_OPD_visit * unit_visit_high) + cost_OPD_tx

keep ScreeningID num* cost*
duplicates drop
save monthly.dta, replace


*** COMBINE BOTH PHASES = OUTPATIENT DATASET ************
use weekly.dta, clear
append using monthly.dta

egen temp_num_OPD_visit = total(num_OPD_visit),  by(ScreeningID)
egen temp_num_OPD_diarr = total(num_OPD_diarr),  by(ScreeningID)

egen temp_cost_OPD_visit = total(cost_OPD_visit), by(ScreeningID)
egen temp_cost_OPD_tx    = total(cost_OPD_tx),    by(ScreeningID)

egen temp_cost_OPD_visit_low = total(cost_OPD_visit_low), by(ScreeningID)
egen temp_cost_OPD_tx_low    = total(cost_OPD_tx_low),    by(ScreeningID)

egen temp_cost_OPD_visit_high = total(cost_OPD_visit_high), by(ScreeningID)
egen temp_cost_OPD_tx_high    = total(cost_OPD_tx_high),    by(ScreeningID)

egen temp_cost_OPD       = total(cost_OPD),       by(ScreeningID)

keep ScreeningID temp*
duplicates drop
rename temp_* *

*** Merge with randomization file to get complete 3300 observations
merge 1:1 ScreeningID using randomization
foreach x of varlist num* cost*{
  replace `x' = 0 if missing(`x') //replace missing values with 0
}
drop _merge

gen yld_diarrhoea = (num_OPD_diarr*dur_diarr_moderate/365.25)*dw_diarr_moderate

save outpatient.dta, replace



***** D. HOSPITAL COSTS + INTERVENTION COSTS = HEALTH SYSTEM COSTS  
*****************************************************************

*** Use outpatient dataset and merge with inpatient dataset
use outpatient.dta, clear
merge 1:1 ScreeningID using inpatient.dta
drop _merge

*** Hospital (health service) costs = Outpatient costs + Inpatient costs
gen cost_hosp  = cost_OPD + cost_IPD

*** Cost of active iron intervention: syrup for arm 1 and mnps for arm 2
generate num_syrp  = (treat == 1)
generate num_mnp    = (treat == 2)
replace  num_mnp    = num_mnp * 90

generate cost_syrp = num_syrp * unit_syrp
generate cost_mnp   = num_mnp   * unit_mnp

*** Cost of programme delivery: for arm 1 and 2 only, not for 3 (placebo)
gen cost_prgm      = 0
replace cost_prgm  = unit_prgm if treat == 1 | treat == 2

*** Cost of intervention (active iron + programme delivery)
gen cost_intv = cost_syrp + cost_mnp + cost_prgm  

*** Total costs (health system) = intervention + hospital
gen cost_hs   = cost_syrp + cost_mnp + cost_prgm + cost_hosp

* create variables for sensitivity analysis
gen cost_hs_syrp_low  = cost_hosp + cost_syrp/1.5 + cost_mnp + cost_prgm
gen cost_hs_syrp_high = cost_hosp + cost_syrp*1.5 + cost_mnp + cost_prgm
gen cost_hs_mnp_low   = cost_hosp + cost_syrp + cost_mnp/1.5 + cost_prgm
gen cost_hs_mnp_high  = cost_hosp + cost_syrp + cost_mnp*1.5 + cost_prgm
gen cost_hs_prgm_low  = cost_hosp + cost_syrp + cost_mnp + cost_prgm/1.5
gen cost_hs_prgm_high = cost_hosp + cost_syrp + cost_mnp + cost_prgm*1.5

gen cost_hs_visit_low  = cost_intv + cost_IPD + cost_OPD_visit_low
gen cost_hs_visit_high = cost_intv + cost_IPD + cost_OPD_visit_high
gen cost_hs_bed_low    = cost_intv + cost_OPD + cost_IPD_bed_low
gen cost_hs_bed_high   = cost_intv + cost_OPD + cost_IPD_bed_high
gen cost_hs_tx_low     = cost_intv + cost_OPD_tx_low  + cost_IPD_tx_low
gen cost_hs_tx_high    = cost_intv + cost_OPD_tx_high + cost_IPD_tx_high


*** Save file and erase files not needed after this step
save healthSystem.dta, replace




***** F. OOP COSTS FROM WEEKLYMORB_HE/MONTHLYMORB_HE DATA
*****************************************************************

*** FOR INTERVENTION PHASE ************
* Use healthcare visits where Diarrhoea was the primary diagnosis
use WeeklyMorb_HE.dta if Diarrhoea == 1, clear

* Using direct non-medical costs to avoid double-counting
keep ScreeningID Q16D Q16E
gen  temp_cost_oop = Q16D + Q16E

* Total costs by child
egen cost_oop      = total(temp_cost_oop), by(ScreeningID)
drop Q16D Q16E temp_cost_oop
duplicates drop
save weekly.dta, replace


*** FOR FOLLOW-UP PHASE ************
* Use healthcare visits where Diarrhoea was the primary diagnosis
use MonthlyMorb_HE.dta if Diarrhoea == 1, clear

* Using only direct non-medical costs to avoid double-counting
// Same reason as above
keep ScreeningID Q10D Q10E
gen  temp_cost_oop = Q10D + Q10E

* Generate total OOP costs by child
egen cost_oop      = total(temp_cost_oop), by(ScreeningID)
replace cost_oop   = cost_oop * 2.17

drop Q10D Q10E temp_cost_oop
duplicates drop
save monthly.dta, replace


* COMBINE BOTH PHASES ************
use weekly.dta, clear
append using monthly.dta

egen temp_cost_oop  = total(cost_oop),  by(ScreeningID)
drop cost_oop
duplicates drop
rename temp_* *

*** Merge with randomization file to get complete 3300 observations
merge 1:1 ScreeningID using randomization.dta
replace cost_oop = 0 if missing(cost_oop)

* Convert from Bangladeshi Taka to US Dollars
replace cost_oop = cost_oop / xch20
drop _merge
sort ScreeningID
save out_of_pocket.dta, replace




***** G. PRODUCTIVITY LOSS FROM WEEKLYMORB_HE/MONTHLYMORB_HE DATA
*****************************************************************

*** FOR INTERVENTION PHASE ************
* Use healthcare visits where Diarrhoea was the primary diagnosis
use WeeklyMorb_HE.dta if Diarrhoea == 1, clear

* Keep only caregivers that reported having to take any time off work/duties
keep if Q17A  == 1

* If only one caregiver took time off work, set Q17B from missing to 0
replace Q17B2 = 0 if missing(Q17B2)

* Add number of hours from all caregivers
gen temp_num_hours = Q17A2 + Q17B2

* Keep only the variables needed
keep ScreeningID temp_num_hours

* Generate total number of hours lost by caregivers of child
egen num_hours  = total(temp_num_hours), by(ScreeningID)
drop temp_num_hours
duplicates drop

*  Save file - productivity loss during the intervention period
save weekly.dta, replace


*** FOR FOLLOW-UP PHASE ************
* Use healthcare visits where Diarrhoea was the primary diagnosis
use MonthlyMorb_HE.dta if Diarrhoea == 1, clear

* Keep only caregivers that reported having to take any time off work/duties
keep if Q11A  == 1

* If only one caregiver took time off work, set Q17B from missing to 0
replace Q11B2 = 0 if missing(Q11B2)

* Add number of hours from all caregivers
gen temp_num_hours = Q11A2 + Q11B2

* Keep only the variables needed
keep ScreeningID temp_num_hours

* Generate total number of hours lost by caregivers of child
egen num_hours  = total(temp_num_hours), by(ScreeningID)
replace temp_num_hours   = temp_num_hours * 2.17
// 2.17 is a factor used to convert 14 days to 1 month. After the intervention perios, children were no longer visited weekly but monthly. During these visits, data were collected based on 14-day recall
drop temp_num_hours
duplicates drop


*  Save file - productivity loss during the follow-up period
save monthly.dta, replace



*** COMBINE BOTH PHASES = OUTPATIENT DATASET ************
use weekly.dta, clear
append using monthly.dta
egen temp_num_hours  = total(num_hours),  by(ScreeningID)
drop num_hours
duplicates drop
rename temp_* *

* Merge with randomization file to get complete 3300 observations
merge 1:1 ScreeningID using randomization.dta
replace num_hours = 0 if missing(num_hours)
drop _merge
sort ScreeningID

* Convert from hours to days
gen num_days = num_hours / 12

* Get GDP per capita per day by dividing GDP per capita by number of days in a year
scalar gdp_pc_pd = 1968.79/365

* Productivity loss = Nuber of days lost * GDP per capita per day
gen cost_prod = num_days * gdp_pc_pd

save productivity.dta, replace



***** H. MERGE COST DATA WITH DEMOGRAPHIC + LABORATORY DATA
*****************************************************************
use healthSystem.dta, clear
merge 1:1 ScreeningID using out_of_pocket.dta
drop _merge
merge 1:1 ScreeningID using productivity.dta
drop _merge

foreach var of varlist cost_hs*{
	gen cost_`var'_ = `var' + cost_oop + cost_prod
	rename *_cost_hs* **
}
rename *_ *


*** Clean up and label variables
format cost* %20.2fc
lab var num_OPD_visit       "# of OPD visits"
lab var num_OPD_diarr       "# of OPD diarrhoea cases"
lab var cost_OPD_visit      "Cost of OPD visits"
lab var cost_OPD_visit_low  "Cost of OPD visits (low value)"
lab var cost_OPD_visit_high "Cost of OPD visits (high value)"
lab var cost_OPD_tx         "Cost of OPD drugs"
lab var cost_OPD_tx_low     "Cost of OPD drugs (low value)"
lab var cost_OPD_tx_high    "Cost of OPD drugs (high value)"
lab var cost_OPD            "Cost of OPD tx"
lab var num_IPD_bed         "# of IPD bed days"
lab var num_IPD_diarr       "# of IPD diarrhoea cases"
lab var cost_IPD_bed        "Cost of IPD bed days"
lab var cost_IPD_bed_low    "Cost of IPD bed days (low value)"
lab var cost_IPD_bed_high   "Cost of IPD bed days (high value)"
lab var cost_IPD_tx         "Cost of IPD drugs"
lab var cost_IPD_tx_low     "Cost of IPD drugs (low value)"
lab var cost_IPD_tx_high    "Cost of IPD drugs (high value)"
lab var cost_IPD            "Cost of IPD tx"
lab var cost_hosp           "Hospital costs"
lab var num_syrp            "# of syrup bottles"
lab var cost_syrp           "Cost of iron syrup"
lab var num_mnp             "# of mnp sachets"
lab var cost_mnp            "Cost of MNPs"
lab var cost_prgm           "Cost of program delivery"
lab var cost_intv           "Cost of intervention"
lab var cost_hs             "Total costs to health system"
lab var cost_hs_syrp_low    "Total costs HS (low syrup value)"
lab var cost_hs_syrp_high   "Total costs HS (high syrup value)"
lab var cost_hs_mnp_low     "Total costs HS (low MNP value)"
lab var cost_hs_mnp_high    "Total costs HS (high MNP value)"
lab var cost_hs_prgm_low    "Total costs HS (low prog delivery value)"
lab var cost_hs_prgm_high   "Total costs HS (high prog delivery value)"
lab var cost_hs_visit_low   "Total costs HS (low visit value)"
lab var cost_hs_visit_high  "Total costs HS (high visit value)"
lab var cost_hs_bed_low     "Total costs HS (low bed value)"
lab var cost_hs_bed_high    "Total costs HS (high bed value)"
lab var cost_oop            "Out of pocket costs"
lab var num_hours           "# of lost productive hours"
lab var num_hours           "# of lost productive days"
lab var cost_prod           "Productivity losses"
lab var cost                "Total costs to society"
lab var cost_syrp_low       "Total costs (low syrup value)"
lab var cost_syrp_high      "Total costs (high syrup value)"
lab var cost_mnp_low        "Total costs (low MNP value)"
lab var cost_mnp_high       "Total costs (high MNP value)"
lab var cost_prgm_low       "Total costs (low prog delivery value)"
lab var cost_prgm_high      "Total costs (high prog delivery value)"
lab var cost_visit_low      "Total costs (low visit value)"
lab var cost_visit_high     "Total costs (high visit value)"
lab var cost_bed_low        "Total costs (low bed value)"
lab var cost_bed_high       "Total costs (high bed value)"

save societal.dta, replace


***** H. MERGE COST DATA WITH DEMOGRAPHIC + LABORATORY DATA
*****************************************************************
 
use Laboratory_HE.dta, clear
rename inwindow14 inwindow14Lab

*** merge with baseline and cost data
merge m:1 ScreeningID using BaseChar_HE.dta
drop  _merge
merge m:1 ScreeningID using societal.dta
drop  _merge

*** drop records laboratory data collected out of window.
*drop if inwindow14Lab == 0
foreach var of varlist HbFgL_venous FerritinngmL fci_total{
    replace `var' = . if inwindow14Lab == 0
}

*** drop variables not needed or are derived variables
drop PID anaemia_mild_venous anaemia_modsev_venous anaemia_mod_venous anaemia_sev_venous fci_total_med IronDA_cappilary inwindow14 inwindow14Lab


*** rename some variables to meaningful names
rename complboth70    compliance
rename Sex_female     female
rename FerritinngmL   ferritin
rename HbFgL_venous   hb
rename anaemia_venous anaemia
rename IronD          ironD
rename IronDA_venous  ironDA

recode     treat 3=0 1=2 2=1
lab define treat 3 "",        modify
lab define treat 0 "Placebo", add
lab define treat 1 "MNP",     modify
lab define treat 2 "Iron",    modify

*** save master dataset
save Analysis_master_data.dta, replace


*** Erase files not needed after this step
cap erase randomization.dta
cap erase weekly.dta
cap erase monthly.dta
cap erase outpatient.dta
cap erase inpatient.dta
cap erase hospital.dta
cap erase healthSystem.dta
cap erase out_of_pocket.dta
cap erase productivity.dta
cap erase societal.dta


