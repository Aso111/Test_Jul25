	*version 12

	* DISCLAIMER: The following code is provided as is and should be used at your own risk. 
	* No user support is available.
	*
	* Nigel Rice, 2011, last updated 30 January 2012
	*
	*
	* This program relates to the material described in Chapters 3 & 12 of the book 
	* Jones, A.M., Rice, N., Bago d'Uva, T. & Balia, S. (2012)
	* "Applied Health Economics. Second Edition", London: Routledge   


	* To run this do-file (main file) in Stata the user must:
	*      1) Install the user-written commands pglm, and fmm using "ssc install" 
	*      and download the command gengam2 available from Anirban Basu's website: 
	*      http://faculty.washington.edu/basua/software.html
	 


	* PRELIMINARIES

	set more off
	clear all
	macro drop _all
	set trace off
	program drop _all

	* Alter the following command to match the folder that you are using:
	cd "c:\My Documents\Routledge book\2nd edition\programs\"

	capture log close
	log using "AHE_2ed_Ch3&12.log", replace

	use "AHE_2ed_Ch_3.dta", clear


	/********************************/
	/* Descriptive analysis of data */

	describe

	summ totexp
	tab posexp
	tab posexp, summarize(posexp)

	drop if totexp == 0
	summ totexp, detail


	/************************************/
	/* PREPARATION OF DATA FOR ANALYSIS */


	clonevar y=totexp
	clonevar lny=ltotexp
	gen sqrty = y^0.5
	global xvars age female income suppins phylim actlim totchr 

	/* Run OLS to remove any item missing data */
	quietly reg totexp $xvars
	gen insamp = 1 if e(sample)  // set sample to estimation sample
	drop if insamp ~= 1

	/* set number of equal groups in data for v-fold cross validation */
	global v=100
	xtile idq = dupersid, nq($v)



	/*************************************************/
	/* Programs to draw upon in evaluation exercise */

	/* Linktest */
	program define lnkt, rclass
	  syntax namelist (min=3) [if] [, ]
	  args cmd y xbhat 
	  tempvar xbhat2
	  gen `xbhat2' = `xbhat'*`xbhat'
	  `cmd' `y' `xbhat' `xbhat2', robust        // general link test
	  test `xbhat2'
	  scalar lnkt = r(p)
	end

	/* Linktest for GLM */
	program define lnktglm, rclass           // Linktest for GLMs
	  syntax namelist (min=3) [if] [, link(string) family(string)]
	  args cmd y xbhat
	  tempvar xbhat2
	  gen `xbhat2' = `xbhat'*`xbhat'
	  `cmd' `y' `xbhat' `xbhat2' , `link' `family' difficult iterate(100)  // general link test
	  test `xbhat2'
	  scalar lnkt = r(p)
	end

	/* Compute R^2, RMSE & MAPE & Pearson's correlation */
	program define evaluate, rclass
	  args a p
	  tempvar pe spe ape yres
	  ereturn clear 
	  qui gen `pe' = (`a'-`p')
	  qui summ `pe'
	  scalar mpe = r(mean)
	  qui gen `spe' = (`a'-`p')^2
	  qui summ `spe'
	  scalar rmse = sqrt(r(mean))
	  qui gen `ape' = abs(`a'-`p')
	  qui summ `ape'
	  scalar mape = r(mean)
	  reg `a' `p'
	  ereturn list
	  scalar r2 = e(r2)
	  gen `yres' = `a'-`p'
	  qui reg `yres' `p'
	  qui test `p'
	  scalar corp = r(p)
	end

	/* Compute modified Hosmer-Lemeshow test  */
	program define hoslem, rclass
	  args a p
	  tempvar yres grp
	  gen `yres' = `a'-`p' 
	  sort `p'
	  xtile `grp' = `p', nq(10)
	  qui tab `grp', gen(gr)
	  qui reg `yres' gr1-gr10, nocons robust
	  qui testparm gr1-gr10
	  scalar hoslm = fprob(e(df_m),e(df_r),e(F))
	  di hoslm
	  drop gr1-gr10
	end

	/* Display in-sample results */
	program define insample, rclass
	 display "In-sample R^2 = " r2
	 display "In-sample RMSE = " rmse 
	 display "In-sample MAPE = " mape
	 display "Pearson correlations (p-value) =" corp
	end

	/* Display Cross-validation results */
	program define crossval
	 display "Cross validation RMSE = " rmse 
	 display "Cross validation MAPE = " mape
	 display "Cross validation MPE = " mpe
	end


	/******************************/
	/* OLS on untransformed costs */

	reg y $xvars, robust
	predict yhat, xb

	/* Linktest */
	qui lnkt reg y yhat
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat 
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample              // calls programme to display results

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui reg y $xvars if idq~=`i' , vce(robust)
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval             // calls programme to display results

	drop yhat yfv


	/******************************/
	/* OLS on transformed costs ***/

	/* Log costs - heteroskedastic errors */

	reg lny $xvars, robust
	linktest
	predict yhat, xb
	gen expr = exp(lny-yhat)
	quietly reg expr yhat
	predict exprhat 
	gen yhatd = exp(yhat) * exprhat 

	/* Linktest */
	qui lnkt reg lny yhat
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhatd  
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhatd
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui reg lny $xvars if idq~=`i' , vce(robust)
	   qui predict yf, xb
	   gen expvf = exp(lny-yf) 
	   quietly reg expvf yf
	   quietly predict expvfhat 
	   gen yfd = exp(yf) * expvfhat 
	   qui replace yfv=yfd if idq==`i'
	   qui drop yf expvf expvfhat yfd  
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yhatd yfv


	/* Square-root costs - heteroskedastic errors */

	reg sqrty $xvars, robust
	predict yhat, xb
	gen rsqr = (sqrty-yhat)^2
	quietly reg rsqr yhat
	predict rsqrhat 
	gen yhatd = (yhat)^2 + rsqrhat 

	/* Linktest */
	qui lnkt reg sqrty yhat
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhatd  
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhatd
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui reg sqrty $xvars if idq~=`i' , vce(robust)
	   qui predict yf, xb
	   gen rsqrvf = (sqrty-yf)^2 
	   quietly reg rsqrvf yf
	   quietly predict rsqrvfhat 
	   gen yfd = (yf)^2 + rsqrvfhat 
	   qui replace yfv=yfd if idq==`i'
	   qui drop yf rsqrvf rsqrvfhat yfd  
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yhatd yfv


	/******************************/
	/* Nonlinear models         ***/

	/**************/
	/* ECM - NLLS */

	gen one = 1
	nl(y=exp({xb: $xvars one})), vce(robust) nolog
	predict yhat if e(sample)

	/* Linktest */
	qui lnkt reg lny yhat
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat 
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui nl(y=exp({xb: $xvars one})) if idq~=`i' , vce(robust)
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf  
	   }
	   
	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yfv


	/********************/
	/* ECM - Poisson-ML */

	poisson y $xvars, vce(robust) nolog
	predict yhat if e(sample)

	/* Linktest */
	qui lnkt reg lny yhat
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat 
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui poisson y $xvars if idq~=`i' , vce(robust)
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf  
	   }
	   
	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yfv


	/*********************/
	/* Generalised Gamma */

	stset y

	/* homoskedastic */

	gengam2 $xvars, vce(robust) nolog time
	predict yxb if e(sample), xb
	predict yhat if e(sample), mean time

	/* Linktest */
	gen yxb2 = yxb^2
	gengam2 yxb yxb2, vce(robust) nolog time
	test yxb2=0
	drop yxb yxb2

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat  
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui gengam2 $xvars if idq~=`i' , vce(robust) nolog time
	   qui predict yf, mean time
	   qui replace yfv=yf if idq==`i'
	   qui drop yf  
	   }
	   
	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yfv

	/* heteroskedastic */

	gengam2 $xvars, ancillary($xvars) vce(robust) nolog time
	predict yxb if e(sample), xb
	predict yhat if e(sample), mean time

	/* Linktest */
	gen yxb2 = yxb^2
	gengam2 yxb yxb2, vce(robust) nolog time
	test yxb2=0
	drop yxb yxb2

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat 
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui gengam2 $xvars if idq~=`i' , ancillary($xvars) vce(robust) nolog time
	   qui predict yf, mean time
	   qui replace yfv=yf if idq==`i'
	   qui drop yf  
	   }
	   
	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yfv

	/*****************************/
	/* Generalised Linear Models */

	/* GLM sqrt-gamma */

	glm y $xvars, l(power 0.5) f(gamma) vce(robust) nolog
	predict yhat if e(sample), mu
	predict yh if e(sample), xb

	/* Linktest */
	linktest, l(power 0.5) f(gamma) vce(robust)
	qui lnktglm glm y yh , link(l(power 0.5)) family(fam(gamma))  
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat 
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui glm y $xvars if idq~=`i', l(power 0.5) f(gamma) vce(robust) nolog
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yh yfv


	/* GLM log-gamma */

	glm y $xvars, l(log) f(gamma) vce(robust) nolog
	predict yhat if e(sample), mu
	predict yh if e(sample), xb

	/* Linktest */
	linktest, l(log) f(gamma) vce(robust)
	qui lnktglm glm y yh , link(l(log)) family(fam(gamma))	
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat  
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui glm y $xvars if idq~=`i', l(log) f(gamma) vce(robust) nolog
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yh yfv


	/* GLM log-normal */

	glm y $xvars, l(log) f(gaussian) vce(robust) nolog
	predict yhat if e(sample), mu
	predict yh if e(sample), xb

	/* Linktest */
	linktest, l(log) f(gamma) vce(robust)
	qui lnktglm glm y yh , link(l(log)) family(fam(gaussian))	
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui glm y $xvars if idq~=`i', l(log) f(gaussian) vce(robust) nolog
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yh yfv

	/* GLM log-poisson */

	glm y $xvars, l(log) f(poisson) vce(robust) nolog
	predict yhat if e(sample), mu
	predict yh if e(sample), xb

	/* Linktest */
	linktest, l(log) f(poisson) vce(robust)
	qui lnktglm glm y yh , link(l(log)) family(fam(poisson))	 
	display "Linktest (p-value) = "lnkt

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat  
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui glm y $xvars if idq~=`i', l(log) f(poisson) vce(robust) nolog
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yh yfv


	/*********************************/
	/* Extended Estimating Equations */

	summ y
	gen scyvar = y/r(mean)
	global sc = r(mean)

	pglm scyvar $xvars
	pglmpredict yhat if e(sample), mu scale($sc)

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat  
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui pglm scyvar $xvars if idq~=`i'
	   qui pglmpredict yf, mu scale($sc)
	   qui replace yfv=yf if idq==`i'
	   qui drop yf
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yfv


	/************************/
	/* Finite Mixture Model */

	fmm y $xvars , vce(robust) components(2) mixtureof(gamma) 
	predict yhat if e(sample)

	/* Modified Hosmer-Lemeshow test */
	qui hoslem y yhat  
	display "Hosmer-Lemeshow test (p-value) = "hoslm

	/* In sample evaluation */
	quietly evaluate y yhat
	insample

	/* v-fold cross validation */
	gen yfv=0
	forvalues i=1/$v {
	   qui fmm y $xvars if idq~=`i', vce(robust) components(2) mixtureof(gamma)
	   qui predict yf
	   qui replace yfv=yf if idq==`i'
	   qui drop yf
	   }

	/* Copas test */
	regress y yfv
	test yfv==1

	quietly evaluate y yfv
	crossval

	drop yhat yfv

	log close 

