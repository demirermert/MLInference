// lassoClean.ado 

capture prog drop lassoClean
program lassoClean, rclass 

syntax varlist(numeric min=1) [if] [in] [, two_way to_indicator(varlist numeric min=1) to_fe(varlist numeric min=1) ]

*(0) Preliminaries	
qui: ssc install isvar
local varlist: list varlist - to_indicator 
local varlist: list varlist - to_fe
local _varlist " "
if wordcount("`two_way'")==1 {
	local max_char 11
}
else {
	local max_char 26
}
foreach var of varlist `varlist' {
	local newname = substr("`var'",1,`max_char')
	qui: gen _`newname'=`var'
	local _varlist "`_varlist' _`newname'"
}

*(1) Turn variables in `to_indicator' into indicators.  
if wordcount("`to_indicator'")!=0 {
local _i " "
foreach var of varlist `to_indicator' {
	qui: replace `var'=. if `var'>.
	local newname = substr("`var'",1,`max_char')
	qui: tabulate `var', missing generate(_`newname'i)
	unab add: _`newname'i*
	local _i "`_i' `add'"
}
}

*(2) Square variables with more than 2 values (any variable with 2 values is collinear with its square)
local _sq " "
foreach var of varlist `_varlist' {
	qui: levelsof `var'
	if wordcount("`r(levels)'")>2 {
		qui: gen `var'_sq= `var'^2
		local _sq "`_sq' `var'_sq"
	}
}

*(3) Missing dummies and missing values
local _mi " "
foreach var of varlist `_varlist' {
	qui: count if missing(`var')
	if `r(N)'>0 {
		qui: gen `var'_mi= (`var'>=.)
		local _mi "`_mi' `var'_mi"
		qui: replace `var'=0 if `var'>=.
	}
}
foreach var of varlist `_sq' {
	qui: replace `var'=0 if `var'>=.
}

*(4) Drop variables that have the same value for everything
local drop " " 
foreach var of varlist _* {
	qui: sum `var'
	if r(min)==r(max) {
		local drop "`drop' `var'" 
		qui: drop `var'
	}
}
local _varlist: list _varlist - drop
local _i: list _i - drop
local _mi: list _mi - drop
local _sq: list _sq - drop

*(5) Get rid of perfectly collinear variables
local drop " "
unab all: _*
local K: word count `all'
forv i=1(1)`=`K'-1' {
	forv j=`=`i'+1'(1)`K' {
		local y: word `i' of `all'
		local x: word `j' of `all'
		qui: isvar `y' `x'
		if wordcount("`r(varlist)'")==2 {
			qui: reg `y' `x'
			if e(r2)==1 {		
				local drop "`drop' `x'" 
				qui: drop `x' 
			}
        }
	}
}
local _varlist: list _varlist - drop
local _i: list _i - drop
local _mi: list _mi - drop
local _sq: list _sq - drop

*(6) Generate two-way interactions
if wordcount("`two_way'")==1 {
local _int " "
unab twoway: `_varlist' `_i' 
local K: word count `twoway'
forv i=1(1)`=`K'-1' {
	forv j=`=`i'+1'(1)`K' {
		local y: word `i' of `twoway'
		local x: word `j' of `twoway'
		local x = substr("`x'",2,.)	
		qui: gen `y'X`x'=`y'*_`x'
		local _int "`_int' `y'X`x'"
	}
}
foreach var of varlist `_int' {
	qui: sum `var'
	if r(min)==r(max) {
		qui: drop `var'
	}
}
}

*(7) Get rid of perfectly collinear variables including interactions
if wordcount("`two_way'")==1 {
unab all: _*
local K: word count `all'
unab non_int: `_varlist' `_i' `_mi' `_sq'
local K1: word count `non_int'
forv i=1(1)`=`K'-1' {
	forv j=`=`i'+1'(1)`K' {
		if `j'>`K1' {
			local y: word `i' of `all'
			local x: word `j' of `all'
			qui: isvar `y' `x'
			if wordcount("`r(varlist)'")==2 {
				qui: reg `y' `x'
				if e(r2)==1 {		
					drop `x' 
				}	
			}
		}
	}
}
}

*(8) Create fixed effects for specified variables 
if wordcount("`to_fe'")!=0 {
foreach var of varlist `to_fe' {
	qui: replace `var'=. if `var'>.
	tabulate `var', missing
	scalar m=length("`r(r)'")
	scalar max=28-m
	local newname = substr("`var'",1,max)
	qui: tabulate `var', missing generate(_`newname'i)
}
}

*(9) Standardize variables
foreach var of varlist _* {
	qui: sum `var'
	qui: replace `var'=(`var'-r(mean))/r(sd)
}

*(10) Output list of lasso-ready controls
unab var_p: _*	
return local var_prepped `var_p'
end




