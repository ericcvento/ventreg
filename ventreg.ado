capture program drop ventreg
program ventreg, byable(recall,noheader)
	syntax [anything] [using/] [if], [, cluster(varlist)]
	args depvar focals rhs
	marksample touse

	*count and process By variables
	if _by() != 0 {	
		local byvars = "`_byvars'"
		local byvarsN : word count("`byvars'")
		
		*FLAG numeric variables****
		forvalues x = 1 / `byvarsN' {
			local byvar`x' : word `x' of `byvars' 
			capture confirm numeric variable `byvar`x'' 
			local byvar`x'isnumeric = _rc
			qui tostring `byvar`x'', replace
		}
			
		preserve
			*Display current by variable values*
			qui keep if `touse'==1
			l `byvars' in 1/1, noobs abbrev(32) 
		restore 	
	}

	*Parse the non-focal indepedent variables into steps
	local steps : word count("`rhs'")
	forvalues x = 1/`steps' {
		local varset`x' : word `x' of `rhs'
	}
	
	*remove the : from variables in the same step
	forvalues y = 1/`steps' {
		tokenize "`varset`y''", parse(:)
		local i = 1
		while "`*'" != "" {
			local var`i' = "`1'"
			macro shift 
			macro shift 
			local i = `i' + 1 
		}
	
		local varset`y' = "`var1'"
		if `i' - 1 > 1 {
			local i = `i' - 1 
			forvalues z = 2/`i' {
				local varset`y' = "`varset`y'' `var`z''" 
			}
		}	
	}
	
	*model one is focal only, load independent variables into local macros
	local model1 = "" 
	forvalues y = 1 /`steps' {
		local yy = `y' + 1 
		local model`yy' "`varset1'" 
		forvalues z = 2/`y' {
			local model`yy' = "`model`yy'' `varset`z''" 
		} 
	}
	
	*CREATE POSTFILE FOR RESULTS*
	if _byindex() == 1 | _byindex() == 0 {
	
		*CREATE DIRECTORY FOR OUTPUT*
		global ventregd = c(current_date)
		global ventregt = c(current_time)
		global ventregt = subinstr("$ventregt",":","_",.)
		global ventregdirname = "C:\Users\Public\Documents\ventreg results_$ventregt"
		mkdir "$ventregdirname"
		sleep 10

		local poststats = "coef se tval pval r2 fval obs focalobs vars focalavg pctfocaldepvar neverfocalavg returncode"
		capture postclose ventregresults

		if _by() == 0 {
			qui postfile ventregresults double(index model) str2000(depvar focalname modelvars) double(`poststats') using "$ventregdirname\ventregresults", replace
		}
		if _by() != 0 {
			if `byvarsN' == 1 {
				qui postfile ventregresults double(index model) str2000(depvar byvars focalname modelvars `byvar1') double(`poststats') using "$ventregdirname\ventregresults", replace
			}
			if `byvarsN' == 2 {
				qui postfile ventregresults double(index model) str2000(depvar byvars focalname modelvars `byvar1' `byvar2') double(`poststats') using "$ventregdirname\ventregresults", replace
			}
			if `byvarsN' == 3 {
				qui postfile ventregresults double(index model) str2000(depvar byvars focalname modelvars `byvar1' `byvar2' `byvar3') double(`poststats') using "$ventregdirname\ventregresults", replace
			}
			if `byvarsN' == 4 {
				qui postfile ventregresults double(index model) str2000(depvar byvars focalname modelvars `byvar1' `byvar2' `byvar3' `byvar4') double(`poststats') using "$ventregdirname\ventregresults", replace
			}
			if `byvarsN' == 5 {
				qui postfile ventregresults double(index model) str2000(depvar byvars focalname modelvars `byvar1' `byvar2' `byvar3' `byvar4' `byvar5') double(`poststats') using "$ventregdirname\ventregresults", replace
			}
			if `byvarsN' == 6 {
				qui postfile ventregresults double(index model) str2000(depvar byvars focalname modelvars `byvar1' `byvar2' `byvar3' `byvar4' `byvar5' `byvar6') double(`poststats') using "$ventregdirname\ventregresults", replace
			}
		}
	}
	
	*REGRESSIONS*
	forvalues y = 1 / `yy' {
		*ENTER LOOP IF AT LEAST ONE OBSERVATION*
		qui count if `touse' == 1 
		if r(N) > 0 {
			assert inlist(`touse',1,0)
			qui count if `touse'==1
			scalar records = r(N)
			scalar regreturncode=999
			if records > 10 {
				*NO CLUSTERING********
				if "`cluster'" == "" {
					capture reg `depvar' `focals' `model`y'' if `touse'
					*qui reg `depvar' `focals' `model`y'' if `touse'==1
				}
				*CLUSTERING*
				else {
					*capture reg `depvar' `focals' `model`y'' if `touse', cluster(`cluster')
					qui reg `depvar' `focals' `model`y'' if `touse'==1, cluster(`cluster')
				}
				scalar regreturncode = _rc
				capture assert inlist(regreturncode,0)
				if _rc!=0 {
					dis "Return Code-" regreturncode}
				}
				
				*warning for missing values by comparing # of records in by group and # of records analyzed*
				qui count if `touse' == 1 
				scalar check1 = r(N) 
				qui count if e(sample) == 1 & `touse' == 1
				scalar check2 = r(N)
				assert records == check1
				if check1 != check2 {
					dis "Missing values for some variables in this model: `model`y''" 
				}
			}

			*****************
			*POST ESTIMATION*
			*****************
			foreach focal in `focals' {
				if regreturncode == 0 {
					scalar coef = _b[`focal']
					scalar stde = _se[`focal']
					scalar tval = coef/stde
					scalar pval = 2*ttail(e(df_r),abs(tval))
					scalar r2 = e(r2)
					scalar fval = e(F)
					scalar prf = Ftail(e(df_m),e(df_r),e(F))
					scalar nvars = e(rank) - 1
					scalar obs = e(N)
					scalar modelvars = e(rank) - 1
					
					*average dependant variable when focal/protected is 1 
					qui sum `depvar' if `focal'==1 & `touse'==1 & e(sample) == 1 
					scalar avgfocaldepvar = r(mean)
					scalar coeffpctdepvar = coef/avgfocaldepvar
					scalar focalobs = r(N)
		
					*AVERAGE DEPVAR IF YOU ARE NEVER A FOCAL*****
					capture confirm variable unprot
					if _rc == 0 {
						dis "CONFLICTING VARIABLE NAME: PLEASE RENAME unprot BEFORE RUNNING THIS PROGRAM"
						exit
					}
					qui gen unprot = 1 if `touse'==1 & e(sample) == 1 
					foreach f in `focals' {
						qui replace unprot = 0 if `f' == 1 
					}
					qui sum `depvar' if unprot==1 & `touse'==1 & e(sample) == 1 
					scalar avgneverfocaldepvar = r(mean)
					drop unprot
				}
				
				else {
					*BLANK ALL STATS (OTHER THAN COUNTS)
					scalar coef = .
					scalar stde = .
					scalar tval = .
					scalar pval = .
					scalar r2 = .
					scalar fval = .
					scalar nvars = .
					scalar modelvars = .
					scalar avgfocaldepvar = .
					scalar coeffpctdepvar = .
					scalar avgneverfocaldepvar = .
					
					qui count if `touse'==1
					assert r(N) == records
					scalar obs = records
					qui count if `focal'==1 & `touse'==1  
					scalar focalobs = r(N)
				}
				
				*LOCAL VARIABLE TO STORE THE NUMERIC VARS*
				local reportstats = "(coef) (stde) (tval) (pval) (r2) (fval) (obs) (focalobs) (modelvars) (avgfocaldepvar) (coeffpctdepvar) (avgneverfocaldepvar) (regreturncode)"
				
				if _by() == 0 {
					post ventregresults (_byindex()) (`y') ("`depvar'") ("`focal'") ("`model`y''") ///
					`reportstats'
				}
				else {
					if `byvarsN' == 1 {
						post ventregresults (_byindex()) (`y') ("`depvar'") ("`byvars'") ("`focal'") ("`model`y''") (`byvar1'[_byn1()]) ///
						`reportstats'
					}
					if `byvarsN' ==  2 {
						post ventregresults (_byindex()) (`y') ("`depvar'") ("`byvars'") ("`focal'") ("`model`y''") (`byvar1'[_byn1()]) (`byvar2'[_byn1()]) ///
						`reportstats'
					}
					if `byvarsN' ==  3 {
						post ventregresults (_byindex()) (`y') ("`depvar'") ("`byvars'") ("`focal'") ("`model`y''") (`byvar1'[_byn1()]) (`byvar2'[_byn1()]) (`byvar3'[_byn1()]) ///
						`reportstats'
					}
					if `byvarsN' ==  4 {
						post ventregresults (_byindex()) (`y') ("`depvar'") ("`byvars'") ("`focal'") ("`model`y''") (`byvar1'[_byn1()]) (`byvar2'[_byn1()]) (`byvar3'[_byn1()]) (`byvar4'[_byn1()]) ///
						`reportstats'
					}
					if `byvarsN' ==  5 {
						post ventregresults (_byindex()) (`y') ("`depvar'") ("`byvars'") ("`focal'") ("`model`y''") (`byvar1'[_byn1()]) (`byvar2'[_byn1()]) (`byvar3'[_byn1()]) (`byvar4'[_byn1()]) (`byvar5'[_byn1()]) ///
						`reportstats'
					}
					if `byvarsN' ==  6 {
						post ventregresults (_byindex()) (`y') ("`depvar'") ("`byvars'") ("`focal'") ("`model`y''") (`byvar1'[_byn1()]) (`byvar2'[_byn1()]) (`byvar3'[_byn1()]) (`byvar4'[_byn1()]) (`byvar5'[_byn1()]) (`byvar6'[_byn1()]) ///
						`reportstats'
					}
				}
			}
		}
	}
	
	*CLOSE POST FILE*
	if _bylastcall() == 1 | _byindex() == 0 {
		postclose ventregresults
		preserve 
			use "$ventregdirname\ventregresults", clear
			qui replace modelvars = trim(modelvars)

			*IF THE BY VARIABLE WAS NUMERIC, RETURN IT.
			if _by() {
				forvalues x = 1 / `byvarsN' {
					if `byvar`x'isnumeric' == 0 {
						qui destring `byvar`x'', replace 
					}
				}
			}
			
			*SAVE TO FILE IN PROGRAM ARGUMENT*
			if "`using'" != "" {
				qui compress
				qui save "`using'", replace
				dis "ventreg, last updated June 5th, 2020."
			}
			else if "`using'"=="" {
			}
			
			*ERASING RESULTS FILE AND DIRECTORY*
			erase "$ventregdirname\ventregresults.dta"
			sleep 10
			rmdir "$ventregdirname"
		restore 
	}
	
	*IF THE BY VARIABLE WAS NUMERIC, RETURN IT.
	forvalues x = 1 / `byvarsN' {
		if `byvar`x'isnumeric' == 0 {
			qui destring `byvar`x'', replace 
		}
	}

end 
