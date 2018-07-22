/********************************************************
PROJECT: RCT etc Paper
PURPOSE: Analyze prepped datasets using LASSO
CODE BY: Maddie McKelway
********************************************************/

clear all
*clear matrix
*clear mata
pause on
set more off
set matsize 11000
set maxvar 15000

local directory "[USER: SET DIRECTORY. SHOULD BE PATH WHERE Lasso FOLDER IS FOUND.]"

cd "`directory'/Analysis/"


*Choose Analysis To Run
local Olken				0
	local Olken0		0
		*These are original tables		
	local Olken1		0
		*This runs lassoClean 
	local Olken2 		1 
		*This runs lasso using pdslasso		

local Duflo				0
	local Duflo0		0
		*These are original tables
	local Duflo1		0
		*This runs lasso clean
	local Duflo2 		1	
		*This runs lasso using pdslasso		
		
	
	
if `Olken'==1 {

if `Olken0'==1 {
global treatment "gen_year1 gen_year1_versi_A"

eststo clear

local outcome1 "pre_natal_visits good_assisted_delivery post_natal_visits iron_pills_categ imm_uptak_pct_23mons times_weighed_last3m vitA_total_6mons_2years mal_weightforage"
foreach x in `outcome1' {
	use "`directory'/Empirical_Examples/Olken_et_al/dta by outcome/`x'", clear
	global control "prov_P_* ind_`x'_w1 ind_`x'_m1 kec_`x'_w1 HH_type_and_panel_* kab_wave_*" 	
	if "`x'"=="pre_natal_visits" | "`x'"== "good_assisted_delivery" | "`x'"== "post_natal_visits" | "`x'"== "iron_pills_categ" | "`x'"== "mal_weightforage"| "`x'"== "enroll_age7to12" | "`x'"== "age7to12_pct_twoweeks" {
		reg `x'_w2 $treatment $control, cluster(Location_ID)
		eststo
	}
	else {
		reg `x'_w2 $treatment $control age_years*, cluster(Location_ID)
		eststo
	}
}
	esttab using "01. Olken/Olken_table3_pt1.tex", keep($treatment) ///
		stats(N, labels("N") fmt(0)) ///
		label nogap b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonotes replace
			
eststo clear	
	
local outcome2 "enroll_age7to12 enroll_age13to15 age7to12_pct_twoweeks age13to15_pct_twoweeks"	
foreach x in `outcome2' {
	use "`directory'/Empirical_Examples/Olken_et_al/dta by outcome/`x'", clear
	global control "prov_P_* ind_`x'_w1 ind_`x'_m1 kec_`x'_w1 HH_type_and_panel_* kab_wave_*" 	
	if "`x'"=="pre_natal_visits" | "`x'"== "good_assisted_delivery" | "`x'"== "post_natal_visits" | "`x'"== "iron_pills_categ" | "`x'"== "mal_weightforage"| "`x'"== "enroll_age7to12" | "`x'"== "age7to12_pct_twoweeks" {
		reg `x'_w2 $treatment $control, cluster(Location_ID)
		eststo
	}
	else {
		reg `x'_w2 $treatment $control age_years*, cluster(Location_ID)
		eststo
	}
}
	esttab using "01. Olken/Olken_table3_pt2.tex", keep($treatment) ///
		stats(N, labels("N") fmt(0)) ///
		label nogap b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonotes replace

}

if `Olken1'==1 {
*First set of outcomes at woman level: pre_natal_visits good_assisted_delivery post_natal_visits iron_pills_categ
use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/all_merged_final.dta", clear
drop if pre_natal_visits==. & good_assisted_delivery==. & post_natal_visits==. & iron_pills_categ==. 
dropmiss _all, force
unab var_raw: ind_pre_natal_visits_w1-k_icp06_i_CDHM
lassoClean `var_raw'
save "`directory'/Empirical_Examples/Olken_et_al/dta merged final/woman_lassoCleaned.dta", replace

*Second set of outcomes at infant level: imm_uptak_pct_23mons times_weighed_last3m vitA_total_6mons_2years mal_weightforage 
use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/all_merged_final.dta", clear
drop if imm_uptak_pct_23mons==. & times_weighed_last3m==. & vitA_total_6mons_2years==. & mal_weightforage==. 
dropmiss _all, force
unab var_raw: ind_pre_natal_visits_w1-k_icp06_i_CDHM
lassoClean `var_raw'
save "`directory'/Empirical_Examples/Olken_et_al/dta merged final/infant_lassoCleaned.dta", replace

*Third set of outcomes at child 7-12 level: enroll_age7to12 age7to12_pct_twoweeks
use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/all_merged_final.dta", clear
drop if enroll_age7to12==. & age7to12_pct_twoweeks==.
dropmiss _all, force
unab var_raw: ind_pre_natal_visits_w1-k_icp06_i_CDHM
lassoClean `var_raw'
save "`directory'/Empirical_Examples/Olken_et_al/dta merged final/child7to12_lassoCleaned.dta", replace

*Fourth set of outcomes at children 13-15 level: enroll_age13to15 age13to15_pct_twoweeks
use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/all_merged_final.dta", clear
drop if enroll_age13to15==. & age13to15_pct_twoweeks==.
dropmiss _all, force
unab var_raw: ind_pre_natal_visits_w1-k_icp06_i_CDHM
lassoClean `var_raw'
save "`directory'/Empirical_Examples/Olken_et_al/dta merged final/child13to15_lassoCleaned.dta", replace

}

if `Olken2'==1 {
	eststo clear
	
	use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/woman_lassoCleaned.dta", clear
	unab controls: _* 	
	global treatment "gen_year1 gen_year1_versi_A"
	
	foreach x in pre_natal_visits good_assisted_delivery post_natal_visits iron_pills_categ {
		pdslasso `x'_w2 $treatment (`controls' prov_P_* HH_type_and_panel_* kab_wave_*), partial(prov_P_* HH_type_and_panel_* kab_wave_*) clus(Location_ID)
			eststo
			estadd local obs_unit "Woman" 
	}
	
	use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/infant_lassoCleaned.dta", clear
	unab controls: _* 
	
	foreach x in imm_uptak_pct_23mons times_weighed_last3m vitA_total_6mons_2years mal_weightforage {
		pdslasso `x'_w2 $treatment (`controls' prov_P_* HH_type_and_panel_* kab_wave_*), partial(prov_P_* HH_type_and_panel_* kab_wave_*) clus(Location_ID)		
			eststo
			estadd local obs_unit "Infant" 			
	}	

	esttab using "01. Olken/Olken_lasso3_pt1.tex", keep($treatment) ///
		stats(obs_unit N, labels("Observation Unit" "N") fmt(0 0)) ///
		label nogap b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonotes replace	
	
	eststo clear 
	
	use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/child7to12_lassoCleaned.dta", clear
	unab controls: _* 	
	
	foreach x in enroll_age7to12 {
		pdslasso `x'_w2 $treatment (`controls' prov_P_* HH_type_and_panel_* kab_wave_*), partial(prov_P_* HH_type_and_panel_* kab_wave_*) clus(Location_ID)
			eststo
			estadd local obs_unit "Child 7-12" 			
	}		
	
	use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/child13to15_lassoCleaned.dta", clear
	unab controls: _* 	

	foreach x in enroll_age13to15 {
		pdslasso `x'_w2 $treatment (`controls' prov_P_* HH_type_and_panel_* kab_wave_*), partial(prov_P_* HH_type_and_panel_* kab_wave_*) clus(Location_ID)
			eststo
			estadd local obs_unit "Child 13-15" 		
	}		
	
	use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/child7to12_lassoCleaned.dta", clear
	unab controls: _* 	
	
	foreach x in age7to12_pct_twoweeks {
		pdslasso `x'_w2 $treatment (`controls' prov_P_* HH_type_and_panel_* kab_wave_*), partial(prov_P_* HH_type_and_panel_* kab_wave_*) clus(Location_ID)
			eststo
			estadd local obs_unit "Child 7-12" 		
	}		
	pause
	use "`directory'/Empirical_Examples/Olken_et_al/dta merged final/child13to15_lassoCleaned.dta", clear
	unab controls: _* 	

	foreach x in age13to15_pct_twoweeks {
		pdslasso `x'_w2 $treatment (`controls' prov_P_* HH_type_and_panel_* kab_wave_*), partial(prov_P_* HH_type_and_panel_* kab_wave_*) clus(Location_ID)		
			eststo
			estadd local obs_unit "Child 13-15" 				
	}
	
	esttab using "01. Olken/Olken_lasso3_pt2.tex", keep($treatment) ///
		stats(obs_unit N, labels("Observation Unit" "N") fmt(0 0)) ///
		label nogap b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonotes replace	
			
}
 
}



if `Duflo'==1 {

if `Duflo0'==1 {
use "`directory'/Empirical_Examples/id_001_educationfertilityhiv/clean_data/studysample_allmerged_lasso.dta", clear
local outcome_2A "dropout05v3 presence evmar05v3 evpreg05v3 ev_p_unm05v3 ev_unp_m05v3"
local outcome_2B "dropout07v2 evmar07v2 evpreg07v2 ev_p_unm07v2 ev_unp_m07v2"
global treatment  "Uonly Honly UH"
global control "yrbirth_all yrbirth_missing date05v3 date07v2 schsize i.stratum"
foreach table in 2A /*2B*/ {
	foreach gender in 2 1 {
		use "`directory'/Empirical_Examples/id_001_educationfertilityhiv/clean_data/studysample_allmerged_lasso.dta", clear
		keep if sex==`gender'
		eststo clear
		foreach outcome in `outcome_`table'' {
			xi: reg `outcome' $treatment $control , cluster(sch03v1)
			eststo 
		}
		esttab using "02. Duflo/Duflo_`table'_`gender'.tex", keep($treatment) ///
			stats(N, labels("N") fmt(0)) ///
			label nogap b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonotes replace			
	}	
	}	
}

if `Duflo1'==1 {
use "`directory'/Empirical_Examples/id_001_educationfertilityhiv/clean_data/studysample_allmerged_lasso.dta", clear

*Drop variables used to form strata
drop zoncode 

*Define lists of controls 
local paper_controls yrbirth_all yrbirth_missing date05v3 date07v2 schsize  
local school_controls date03v1 kcpe1997 kcpe2000 kcpe2001 kcpe2002 kcpe2003 kcpe2004 boys2002 girls2002 total02 ratio02 sdkcpe ///
	Nnew femaleht Necd Nlower Nupper Nuplo Npta Nfemale Nfemaleprim ///
	TOTteachers meanage meanexp latrine_2004 situation ///
	total_2km total_3km total_4km total_5km total_6km total_8km	
local to_dummy discode divcode pres03v1 
	/*Notes: 
		(1) training variables (trainingsession Nfemtraining Nmaletraining Nteacherstraining trainingstart_date Ntrained Nfemtrained Nupperfemaletrained Nuplofemaletrained)
			describe training teachers got in the gov't curriculum. These variables are almost perfectly correlated with treatment 
			assignment so shouldn't be potential controls. In addition, they aren't exactly baseline covariates - they same something
			about how treatment was implemented. 
		(2) htpresent dhtpresent not observed for control or Uonly. 	
		*/
		
*Run lasso clean
drop _*
lassoClean `paper_controls' `school_controls', two_way to_indicator(`to_dummy')
save "`directory'/Empirical_Examples/id_001_educationfertilityhiv/clean_data/studysample_allmerged_lassoCleaned.dta", replace

}

if `Duflo2'==1 {
use "`directory'/Empirical_Examples/id_001_educationfertilityhiv/clean_data/studysample_allmerged_lassoCleaned.dta", clear
qui: tabulate stratum, generate(stratum_i)
global strata "stratum_i*"
unab prepped: _* 

local outcome_2A "dropout05v3 presence evmar05v3 evpreg05v3 ev_p_unm05v3 ev_unp_m05v3"
local outcome_2B "dropout07v2 evmar07v2 evpreg07v2 ev_p_unm07v2 ev_unp_m07v2"
global treatment  "Uonly Honly UH"

foreach table in 2A /*2B*/ {
	foreach gender in 2 1 {
		preserve 
		keep if sex==`gender'
		eststo clear
		foreach outcome in `outcome_`table'' {		
			pdslasso `outcome' $treatment (`prepped' $strata), partial($strata) cluster(sch03v1)
				eststo	
				if "`gender'"=="1" {
					estadd local sample "Boys"
				}
				if "`gender'"=="2" {
					estadd local sample "Girls"
				}				
		}
		if "`table'"=="2A" {		
			esttab using "02. Duflo/Duflo_LASSO_`table'_`gender'.tex", keep($treatment) ///
				stats(sample N, labels("Sample" "N") fmt(0 0)) ///
				label nogap b(3) se(3) star(* 0.10 ** 0.05 *** 0.01) nonotes replace				
		}
		restore 
	}
}	
}
}




