program define ppml_fe_bias, eclass

*! Stata package for correcting inferences from two-way and three-way FE-PPML models
*! contact: Tom Zylkin
*! Department of Economics, University of Richmond
*! This version: v1.1, June 2020
*!
*! Suggested citation: Weidner, Martin and Thomas Zylkin (2020): 
*! "Bias and Consistency in Three-way Gravity Models"
*! arXiv preprint arXiv:1909.01327.

version 13.1
syntax varlist [if] [in], ///
		lambda(varname)  			///
		i(varname)                  ///
		j(varname)                  ///
	[	                            ///
	    t(varname)                  ///
		bias(name)                  ///
		v(name)						///
		w(name)					    ///
		b_term(name)				///
		d_term(name)                ///
		b1(string)   				///
		b2(name)					///
		d1(name)					///
		d2(name)					///
		beta(name)					///
		NOSTERR						///
		approx                      ///
		exact						///
		twoway                      ///
		NOTABLE						///
	]


// This program makes use of four excellent user-created packages:
//   - hdfe by Sergio Correia
//   - gtools by Mauricio Caceres Bravo
//   - rowmatutils by Matthew Baker
//   - frmttable by John Luke Gallup

cap which frmttable
if _rc == 111 {
	di in red "You will need to install -outreg- in order to use this command."
	di
	di in red "To install, type -ssc install outreg-".
	exit 111

}

cap which hdfe
if _rc == 111 {
	di in red "You will need to install -hdfe- in order to use this command."
	di
	di in red "To install, type -ssc install hdfe-".
	exit 111
}

cap which gtools
if _rc == 111 {
	di in red "You will need to install -gtools- in order to use this command."
	di
	di in red "To install, type -ssc install gtools-".
	exit 111
}

cap findfile lrowmat_utils.mlib
if _rc == 601 {
	di in red "You will need to install -rowmat_utils- in order to use this command."
	di
	di in red "To install, type -ssc install rowmat_utils-".
	exit 111
}

if "`twoway'" == "" & "`t'"=="" {
	di in red "A name for the t() variable must be provided unless the -twoway- option is also enabled"
}
if "`twoway'" != "" & "`t'"=="" {
	tempvar t
	gen `t' =1
	di in red "Since no time id has been provided, it has been assumed that there only is one time period."
}

tokenize `varlist'
local y `1'
macro shift
local xvars `*'
unab  xvars: `*'

tempvar sum_lambda pair_id t_id

tempname demeaned

local names "w b_term d_term b1 b2 d1 d2"
foreach n of local names {	
	if "``n''" == "" {
		tempname `n'
	}
}

marksample todo

qui hashsort `i' `j' `t'
by `i' `j': gegen `sum_lambda'_ij = sum(`lambda') if `todo'

qui replace `todo' = 0 if `sum_lambda'_ij == 0
qui replace `todo' = 0 if `lambda' == 0 

gegen `pair_id' = group(`i' `j') if `todo'

gegen `t_id' = group(`t')        if `todo'     // time_id must be 1,2,3...

// make sure all id's are numeric.
capture confirm numeric `i'
if !_rc {
	tempvar i_var
	gegen `i_id' = group(i) if `todo'
}
else{
	local i_id = "`i'"
}
capture confirm numeric `j'
if !_rc {
	tempvar j_var
	gegen `j_id' = group(j) if `todo'
}
else {
	local  j_id  = "`j'"
}

if "`exact'" == "" & "`approx'" == "" & "`twoway'" == ""  {
	
	qui gdistinct `i_id'
	local n_i = r(ndistinct)
	
	qui gdistinct `j_id'
	local n_j = r(ndistinct)
	
	qui gstats sum `t_id'
	local n_t = r(max)
	
	local n_twoway_fes = `n_i'*`n_t' + `n_j'*`n_t'
	
	if `n_twoway_fes' > 1500 {
		local approx = "approx"
	}
}

if "`twoway'" != ""{
	local approx = ""
	local nosterr = ""
}

if "`exact'" != ""{
	local approx = ""
}


// only need this objects if using an approximation for the standard error corrections.
if "`approx'" != "" {
	sort `i' `t' `j' 
	by `i' `t': gegen `sum_lambda'_it = sum(`lambda') if `todo'
	
	sort `j' `t' `i' 
	by `j' `t': gegen `sum_lambda'_jt = sum(`lambda') if `todo'
	
	di in red "note: because of the size of the data, an approximation will be used to compute the adjusted variance. Use the -exact- option if you wish to compute the variance exactly."
}


// weighted demeaning of X wrt FEs, weighted by lambda
if "`twoway'" != "" {
	qui hdfe `xvars' if `todo' [aw=`lambda'], absorb(`i_id'#`t_id' `j_id'#`t_id') gen(`demeaned')
}
else {
	qui hdfe `xvars' if `todo' [aw=`lambda'], absorb(`i_id'#`t_id' `j_id'#`t_id' `i_id'#`j_id') gen(`demeaned')
}

local before : list sizeof xvars
qui _rmcoll `demeaned'* if `todo' & `lambda'>0, forcedrop  // check simple collinearity across `policyvars'
local okvars = r(varlist)
local after : list sizeof okvars
if `before' != `after' {
	di in red "The set of x variables (`xvars') does not appear to be of full rank after conditioning on the fixed effects."
	exit 111
}

if "`bias'" == "" {
	tempname bias
}
if "`v'" == "" {
	tempname v
}
tempname orig_ses

if "`twoway'" == "" {
	mata: analyticalbiascorrection("`demeaned'*", "`i_id'", "`j_id'", "`t_id'", "`pair_id'", "`y'", "`lambda'", "`sum_lambda'", "`bias'", "`v'", ///
									"`w'", "`b_term'", "`d_term'", "`b1'", "`b2'", "`d1'", "`d2'", "`approx'", "`orig_ses'", "`nosterr'",  "`todo'" )  
}
else {
	mata: twowayse_correction("`demeaned'*", "`i_id'", "`j_id'", "`t_id'", "`pair_id'", "`y'", "`lambda'", "`v'", "`w'", "`orig_ses'", "`todo'" )
}

// what type of table to make?
if "`notable'" != "" {
	local table =0
}
else if("`nosterr'"=="" & "`beta'" !="") {
	local table = 1
}
else if("`nosterr'"=="" & "`beta'" =="") {
	local table = 2
}
else if ("`nosterr'"!="" & "`beta'" !="") {
	local table = 3
}
else if ("`nosterr'"!="" & "`beta'" =="") {
	local table = 4
}

local n_vars : word count `xvars'

// make sure that user-input beta is configured as a row vector and has the right number of elements.
if "`beta'" != "" & "`notable'" == "" {
	if colsof(`beta') > rowsof(`beta') {
		matrix `beta' = `beta''
	}
	if rowsof(`beta') > `n_vars' {
		matrix `beta' = `beta'[1..`n_vars',1]
		local beta_name = "`beta'"
		di in red "note: `beta_name' matrix will be shortened to the same length as the number of x-variables"
	}
	else if rowsof(`beta') < `n_vars' {
		local table = `table' + 1
	}
}	

qui sum `t_id'
if `r(max)'==1 {
	local se_note = "Robust standard errors"
}
else {
	local se_note = "Standard errors clustered by pair"
}

tempname results_matrix
if `table' == 1 {
	ereturn clear
	tempname ses adj_betas store_results
	mata  st_matrix("`ses'", sqrt(diagonal( st_matrix("`v'") )))
	
	if ("`twoway'" != "") {
		matrix `adj_betas' = `beta'
	}
	else{
		matrix `adj_betas' = `beta' - `bias'
	}
	matrix rownames `adj_betas' = `xvars'
	matrix rownames `v' = `xvars' 
	matrix colnames `v' = `xvars'
			
	matrix  `results_matrix' = `adj_betas', `ses'
	
	matrix rownames `results_matrix' = `xvars'
	if ("`twoway'" == "") {
		matrix rownames `bias' = `xvars'
	}
	matrix rownames `beta' = `xvars'
	matrix rownames `ses'  = `xvars'
	matrix rownames `orig_ses' = `xvars'
	
	tempname orig_results
	matrix `orig_results' = `beta', `orig_ses'
	
	qui frmttable, statmat(`orig_results') substat(1) sdec(7) ctitles("", "original")  
	if ("`twoway'" == "") {
		qui frmttable, statmat(`bias') sdec(7) ctitles("", "bias")  merge
	}
	qui frmttable, statmat(`ses')  sdec(7) ctitles("", "adjusted SEs")  merge
	
	local bc = rowsof(`results_matrix')
    matrix stars = J(`bc',2,0)
    forvalues k = 1/`bc' {
        matrix stars[`k',2] =   ///
       (abs(`results_matrix'[`k',1]/`results_matrix'[`k',2]) > abs(invnormal(0.10/2))) +   ///
	   (abs(`results_matrix'[`k',1]/`results_matrix'[`k',2]) > abs(invnormal(0.05/2))) +   ///
       (abs(`results_matrix'[`k',1]/`results_matrix'[`k',2]) > abs(invnormal(0.01/2)))
	}
	frmttable, statmat(`results_matrix') substat(1) sdec(7) ctitles("", "bias-corrected") ///
	           note("`se_note', using a local de-biasing adjustment"\ ///
			        " to account for estimation noise in the `i'-`t' and `j'-`t' fixed effects."\ ///
					"* p<0.10; ** p<0.05; *** p<0.01") ///
			   annotate(stars) asymbol(*,**,***)  merge
	
	matrix `adj_betas' = `adj_betas''
	matrix `v' = (`v' + `v'')/2
	ereturn post `adj_betas' `v', depname(`y')
	ereturn local cmdline "ppml_fe_bias `0'"
	ereturn local cmd "ppml_fe_bias"
}
if `table' == 2 {
	tempname ses
	mata  st_matrix("`ses'", sqrt(diagonal( st_matrix("`v'") )))

	if ("`twoway'" == "") {
		matrix `results_matrix' = `bias', `ses'
		matrix rownames `results_matrix' = `xvars'
		frmttable, statmat(`results_matrix') sdec(7) ctitles("", "bias", "adjusted SE")  ///
					note("`se_note', using a local de-biasing adjustment"\ ///
			        "to account for estimation noise in the `i'-`t' and `j'-`t' fixed effects.")
	}
	else{
		matrix rownames `ses' = `xvars'
		frmttable, statmat(`ses') sdec(7) ctitles("", "adjusted SE")  ///
					note("`se_note', using a local de-biasing adjustment"\ ///
			        "to account for estimation noise in the `i'-`t' and `j'-`t' fixed effects.")
	}
}
if `table' == 3 {
	tempname adj_betas
	if ("`twoway'" == "") {
		matrix `adj_betas' = `beta'
	}
	else{
		matrix `adj_betas' = `beta' - `bias'
	}
	matrix `results_matrix' = `beta'', `bias', `adj_betas' 
	matrix rownames `results_matrix' = "`xvars'"
	
	frmttable, statmat(`results_matrix') sdec(7) ctitles("", "original", "bias", "bias-corrected")
}
if `table' == 4 {
	matrix `results_matrix' = `bias'
	matrix rownames `results_matrix' = "`xvars'"
	frmttable, statmat(`results_matrix') sdec(7) ctitles("", "bias")
}							   
end

mata: 
scalar altsum (vector r, real todo) {
	return(sum(r:*todo))
}
end


mata: 
void analyticalbiascorrection(string scalar demeaned_x, string scalar i_var, 
				              string scalar j_var, string scalar t_var, string scalar pair_var, string scalar y_var, string scalar lam, 
				              string scalar sum_lam, string scalar bias_name, string scalar V_name, string scalar W_name, string scalar B_name, string scalar D_name,
							  string scalar B1_name, string scalar B2_name, string scalar D1_name, string scalar D2_name, string scalar approx_lev, string scalar orig_sterr, string scalar no_sterr,
							  | string scalar todo_var)
{
	if(approx_lev=="approx") {
		vars = st_data(., (demeaned_x), todo_var)
		vars = st_data(., (i_var, j_var, t_var, pair_var, y_var, lam, sum_lam+"_ij",sum_lam+"_it",sum_lam+"_jt"), todo_var), vars
		K   = cols(vars)-9
	}
	else {
		vars = st_data(., (demeaned_x), todo_var)
		vars = st_data(., (i_var, j_var, t_var, pair_var, y_var, lam, sum_lam+"_ij"), todo_var), vars
		K   = cols(vars)-7
	}
	
	T         = max(vars[.,3])
	NN_panels = max(vars[.,4]) 
	NNT = NN_panels * T
	
	ij = vars[.,1] + vars[.,2]/(2*max(vars[.,2]))
	ij = uniqrows(ij)#J(T,1,1)
	i  = floor(ij)
	j  = floor((ij-i):*(2*max(vars[.,2])) :+ 1e-6)
	
	t  = J(NN_panels,1,1)#(1..T)'
	it = i + t/(2*T)
	jt = j + t/(2*T)

	uniq_i = uniqrows(i)
	uniq_j = uniqrows(j)
	uniq_t = (1..T)'
	uniq_it = uniqrows(it)
	uniq_jt = uniqrows(jt)
	
	N_i = rows(uniq_i)
	N_j = rows(uniq_j)
	
	N_j_1 = NN_panels/N_i
	N_i_1 = NN_panels/N_j
	
	if(approx_lev=="approx") {
		// fill out panels with missing years
		index =  T :* (vars[.,4]:-1) :+ vars[.,3]
		tempvars = J(NNT, 5+K, 0)
		tempvars[index,.] = editmissing(vars[., 5..(9+K)],0)
	}
	else {
		// fill out panels with missing years
		index =  T :* (vars[.,4]:-1) :+ vars[.,3]
		tempvars = J(NNT, 3+K, 0)
		tempvars[index,.] = editmissing(vars[., 5..(7+K)],0)
	}
	
	y             = tempvars[.,1]
	lambda        = tempvars[.,2]
	sum_lambda_ij = tempvars[.,3]
	if(approx_lev=="approx") {
		sum_lambda_it = tempvars[.,4]
		sum_lambda_jt = tempvars[.,5]
		demeanedX     = tempvars[.,6..(5+K)]
		theta_it = editmissing(lambda:/sum_lambda_it,0)
		theta_jt = editmissing(lambda:/sum_lambda_jt,0)
		theta_ij = editmissing(lambda:/sum_lambda_ij,0)
	}
	else {
		demeanedX     = tempvars[.,4..(3+K)]
		theta_ij = editmissing(lambda:/sum_lambda_ij,0)
	}
	
	e = y - lambda
	
	W = (1/NN_panels) * ((lambda:*demeanedX)' * demeanedX) 
	
	
	// have to assume all t are populated (balanced panels) for this to work
	s_equals_t = J(NN_panels,1,1)#I(T)
	
	// t_equals_s is an (NNT) x T matrix with entries equal to 1 when s=t, where s indexes columns.
	
	temp1 = theta_ij#J(1,T,1)
	temp2 = colshape(theta_ij, T) #J(T,1,1)
	temp = temp1:*temp2	
	H_ij = (-temp :+ theta_ij:*s_equals_t):*sum_lambda_ij  

	// H_ij is a matrix that stacks ij-specific TxT blocks, with each ij-specific block containing H_ij.
	
	/*
	// G: Would be stored in memory as an (NNTT) x (T) matrix. Each (TT)xT block would be a single ij-specific element of G 
	
	// each TxT block will be indexed by r
	theta_r = theta_ij # J(T,1,1)
	
	// within each TxT block in G, rows will be indexed by s and columns will be indexed by t
	theta_s = colshape(colshape(theta_ij,T)#J(T,1,1),1)#J(1,T,1)
	
	theta_t = colshape(theta_ij,T)#J(T*T,1,1)
	
	r_equals_t = (t # J(T, T, 1) :== uniq_t')
	
	s_equals_t = (colshape(colshape(t,T)#J(T,1,1),1) # J(1,T,1) :== uniq_t')
	
	r_equals_s = (t # J(T, T, 1)) :== (colshape(colshape(t,T)#J(T,1,1),1))


	G = /*-2 :* theta_r :* theta_s :* theta_t     + */
	          r_equals_s :* theta_r :* theta_t  +
			  r_equals_t :* theta_r :* theta_s  +
			  s_equals_t :* theta_s :* theta_r  /*-
			  r_equals_t :* s_equals_t :* theta_r */
		
	G = G :* (sum_lambda_ij#J(T,1,1))
	*/

	// Fast trace is just eg if T=3 column 1 + column 5 + column 9
	fast_trace = rowshape(I(T), T*T)

	i_short = colshape(i,T)[.,1]
	j_short = colshape(j,T)[.,1]
	
	collapse_i = ((i_short#J(1,N_i,1))':==uniq_i)

	collapse_j = ((j_short#J(1,N_j,1))':==uniq_j)
	
	/* create separate Mata function with fast method; link using
	   if statement. */

	if (no_sterr == "") {	
		if (approx_lev=="approx") {
			
			temp_it = colshape(theta_it,T)#J(T,1,1) :* J(NN_panels,1,1)#I(T)		 
			temp_jt = colshape(theta_jt,T)#J(T,1,1) :* J(NN_panels,1,1)#I(T)
			temp_ij = colshape(theta_ij,T)#J(T,1,1) 
			
			/*
			DDDD = -WM_it + WM_it*WM_ij - WM_jt + WM_jt*WM_it - WM_jt*WM_it*WM_ij + 
					WM_jt*WM_ij - WM_ij + WM_ij*WM_it - WM_ij*WM_it*WM_ij + WM_ij*WM_jt - 
					WM_ij*WM_jt*WM_it + WM_ij*WM_jt*WM_it*WM_ij - WM_ij*WM_jt*WM_ij */
		
			// 1, 3, 4, 7, 8, 10, 11
			temp_DDDD = -temp_it - temp_jt + temp_jt:*temp_it - temp_ij + colshape(theta_ij :*theta_it,T)#J(T,1,1) +
						colshape(theta_ij:*theta_jt,T)#J(T,1,1) - 
						colshape(theta_ij:*rowsum(temp_jt:*temp_it),T)#J(T,1,1)

			stacked_eye = J(NN_panels,1,1)#I(T)
			opp_eye     = J(NNT,T,1)-stacked_eye
			temp_DDDD = temp_DDDD + colshape(theta_it:*theta_ij,T)#J(T,1,1) :* stacked_eye
			
			for (k=1; k <= T-1; k++) {	
				temp_DDDD = temp_DDDD + (theta_it[1..NNT-k]:*theta_ij[1+k..NNT] \ J(k,1,0)) :* (J(NNT,k,0),stacked_eye[.,1..T-k]) + 	
										(J(k,1,0)\theta_it[1+k..NNT]:*theta_ij[1..NNT-k] )  :* (stacked_eye[.,1+k..T],J(NNT,k,0))
			}
	
			temp5 = temp_jt:*temp_it:*theta_ij
			temp6 = temp_jt:*theta_ij

			for (k=1; k <= T; k++) {
				temp5[.,k] = temp5[.,k]+theta_jt:*theta_it:*temp_ij[.,k] :* opp_eye[.,k]
				temp6[.,k] = temp6[.,k]+theta_jt:*temp_ij[.,k] :* opp_eye[.,k]
			}
			temp_DDDD = temp_DDDD - temp5 + temp6
			
			index = (1..NN_panels)*T

			for (k=1; k <= T; k++) {
				temp_DDDD[.,k] = temp_DDDD[.,k] - (rowsum(colshape(theta_ij :* theta_it,T)) :* theta_ij[index:-(T-k),.])#J(T,1,1)
			}
			for (k=1; k <= T; k++) {
				temp_DDDD[.,k] = temp_DDDD[.,k] + rowsum(colshape(theta_ij :* temp5[,k],T))#J(T,1,1)
			}
			for (k=1; k <= T; k++) {
				temp_DDDD[.,k] = temp_DDDD[.,k] - (rowsum(colshape(theta_ij :* theta_jt,T)) :* theta_ij[index:-(T-k),.])#J(T,1,1)
			}
		}	
		else{

			// sort of a fast xi / egen to construct within-transformed sets of it- and jt- dummies
			temp = s_equals_t - colshape(theta_ij,T)#J(T,1,1)

			/*
			d_t   = (t # J(1, T, 1) :== uniq_t')
			temp = colshape(rowsum(theta_ij :* d_t),T)#J(T,1,1) 
			temp = d_t - temp
			*/
			
			d_it_tilde = J(NNT,N_i*T,0)
			d_jt_tilde = J(NNT,N_j*T,0)

			for (k=1; k <= N_i; k++) {
				idx1 = (1, (k-1)*T+1 \ NNT ,k*T)
				idx2 = (i:==uniq_i[k]) 
				d_it_tilde[|idx1|] = temp:*idx2
			}

			for (k=1; k <= N_j; k++) {
				idx1 = (1, (k-1)*T+1 \ NNT ,k*T)
				idx2 = (j:==uniq_j[k]) 
				d_jt_tilde[|idx1|] = temp:*idx2
			}	

			d_tilde = (d_it_tilde, d_jt_tilde)
			d_ijt = d_tilde :> 0
	
			// This can be shown to be equal to inv (sum_ij (d_ij' H_ij d_ij))
			V_FE = invsym((d_tilde:*lambda)' * (d_tilde))

			// Next, we need d_ij V_FE d_ij' (though actually we only need the TxT blocks that lie along the diagonal!)
			//   - d_ijt  is NNT x F  (NN vertically stacked T x F matrices)
			//   - V_FE   is F x F    
			//   - d_ijt' is F x NNT  (NN horizontally stacked F x T matrices)

			// d_ijt x V_F gives me another NNT x F matrix. 
			// - the elements are sum_f {d_ij1,f V_f1} sum_f {d_ij1,f V_f2} sum_f {d_ij1,f V_f3} ...
			//                    sum_f {d_ij2,f V_f1} sum_f {d_ij2,f V_f2} sum_f {d_ij2,f V_f3} ... 
			//                    sum_f {d_ij3,f V_f1} sum_f {d_ij3,f V_f2} sum_f {d_ij3,f V_f3} ... 
	
			// Now consider (d_ijt x V_F) :* d_ijt
			// - the elements would be sum_f {d_ij1,f V_f1} d_ij1,1 sum_f {d_ij1,f V_f2}d_ij1,2 sum_f {d_ij1,f V_f3} d_ij1,3 ...
			//                         sum_f {d_ij2,f V_f1} d_ij2,1 sum_f {d_ij2,f V_f2)d_ij2,2 sum_f {d_ij2,f V_f3} d_ij2,3 ... 
			//                         sum_f {d_ij3,f V_f1} d_ij3,1 sum_f {d_ij3,f V_f2)d_ij3,2 sum_f {d_ij3,f V_f3} d_ij3,3 ... 
	
			// Finding the rowsum then gives you the main diagonal of d_ij V_FE d_ij'
			// To get the full outer product, pull apart the lefthand side so that it repeats.
			// For the right-hand side, reshape, pull apart, then reshape again.
	
			dVd  = J(NN_panels, T*T,0)
			dV    = d_ijt * V_FE
			dVd[.,(0..(T-1))*T :+ (1..(T))] = colshape(rowsum(dV :* d_ijt),T)
	
			for (u=2; u <= T; u++) {
				dVd[.,(0..(T-u))*T :+ (u..(T))] = colshape(rowsum(dV[selectindex(t:<(T-u+2)),.] :* d_ijt[selectindex(t:>(u-1)),.]),T-(u-1))	
		
				if(u==2) {
					idx1 = ((u-1)..(T-1))*T:+(u-1)
					idx2 = (u-2):+(u..T)
				}
				else{
				idx1 = idx1,((u-1)..(T-1))*T:+(u-1)
				idx2 = idx2,(u-2)*T:+(u..T)
				}
			}
			dVd[.,idx1] = dVd[.,idx2]  // this is the NN TxT blocks that lie aloneg the diagonal of d_ij V_FE d_ij'; each TxT block is laid out as a (TxT)x1 vector.
		}
	}

	// This gets us x_ij V_x x_ij'
	V_X = (1/NN_panels) * invsym(W)
	
	tempX = colshape(demeanedX, T*K)#J(T,1,1)
	tempX = colshape(tempX, K) 
	
	xVx = rowsum( ((demeanedX * V_X) # J(T,1,1)) :* tempX)
	
	xVx = colshape(xVx, T*T)
	
	// sum_j H_ij and sum_i H_ij (then invert and compute trace)
	H_i = collapse_i * colshape(H_ij,T*T)  

	H_j = collapse_j * colshape(H_ij,T*T)  
	
	//Fast inversion can be accomplished via rm_newtinv() 
	H_i_inv = rm_newtinv(H_i,30,1e-12)   // can confirm using: pinv(colshape(H_i,T)[|1,1 \T,T|])   /* uses Moore-Penrose pseudoinverse */
	
	H_j_inv = rm_newtinv(H_j,30,1e-12)
	
	// SS_ij is the OUTER product of scores
	SS_ij = (e # J(1,T,1)) :*(colshape(e,T)#J(T,1,1))
	SS_ij = colshape(SS_ij, T*T)

	SS_i = (collapse_i * SS_ij)

	SS_j = (collapse_j * SS_ij)
	
	if (no_sterr == "") {
		if (approx_lev=="approx") {
			lev_correction = -rm_transpose(colshape(temp_DDDD,T*T)) + rm_matmult(colshape(H_ij,T*T), xVx)
			lev_correction = rm_newtinv(colshape(fast_trace,T*T):-lev_correction,30,1e-12)
		}
		else{
			lev_correction = rm_matmult(colshape(H_ij,T*T), xVx+dVd)
			lev_correction = rm_newtinv(colshape(fast_trace,T*T):-lev_correction,30,1e-12)
		}
		SSh_ij = rm_matmult(lev_correction, SS_ij)
		
		/*
		SSh_i = (collapse_i * SSh_ij)
		SSh_j = (collapse_j * SSh_ij)
		*/
	}
	
	// construct OMEGA^U
	// - demeanedX  is NNT x K (think of as NN T x K matrices, stacked vertically)
	// - demeanedX' is K x NNT (think of a NN K x T matrices, stacked horizontally)
	// - SSh_ij is NN x TT (think of as NN TT row vectors, stacked vertically)
	
	tempX     = colshape (demeanedX, K*T)
		
	X_reshape = J(1,K,1)#(1..T) :* K :+ (-(K-1)..0)#J(1,T,1)

	tempX = colshape(tempX[.,X_reshape],T)           // reshapes demeanedX so that it is now NNK x T - or NN K x T matrices vertically stacked.	
	
	if (no_sterr == "") {	
		tempO = rm_matvecmult(SSh_ij#J(K,1,1), tempX)    // gives me NN K x T matrices, stacked vertically
		
		// next I want to construct NN KxK outer products 
		tempO = (tempX#J(K,1,1)) :* colshape(colshape(tempO,T*K)#J(K,1,1),T) 
			// (I)  NN blocks; within each block, there are KK rows of T columns, with each k repeated K times
			// (II) NN blocks; within each block, there are KK rows of T columns, with 1..K repeated K times 
	
		tempO = rowsum(tempO)
	
		tempO = colsum(colshape(tempO, K*K))
		
		OMEGA = (1/(NN_panels)) * colshape(tempO,K)
	
		V     = (1/(NN_panels)) * invsym(W) * OMEGA * invsym(W)
	
		//  Stata uses a finite-sample correction to calculate SEs. See: 
		// 	- https://www.stata.com/meeting/13uk/nichols_crse.pdf.
		//  - https://www.stata.com/manuals/u20.pdf (p. 52)
		//	- https://www.stata.com/manuals/p_robust.pdf (p. 13)
		V = NN_panels / (NN_panels-1) * V
		
		"Adjusted SEs"
		sqrt(diagonal(V))
		
		st_matrix(V_name, V)
	}

	//original VCV matrix
	tempO = rm_matvecmult(SS_ij#J(K,1,1), tempX)
	tempO = (tempX#J(K,1,1)) :* colshape(colshape(tempO,T*K)#J(K,1,1),T)
	tempO = rowsum(tempO)
	tempO = colsum(colshape(tempO, K*K))
	OMEGA = (1/(NN_panels)) * colshape(tempO,K)
	V_orig = (1/(NN_panels)) * invsym(W) * OMEGA * invsym(W)
	V_orig = NN_panels / (NN_panels-1) * V_orig
	
	//original SEs (for comparison)
	st_matrix(orig_sterr, (diagonal(sqrt(V_orig))))

	B = D = J(K,1,0)
	B1 = B2 = J(N_i,K,0)
	D1 = D2 = J(N_j,K,0)
	
	tempH = colshape(H_ij,T*T)
	//tempG = colshape(G,T*T)

	for (k=1; k <= K; k++) {
		
		// compute xHS objects that appear in the first terms in the bias.
		tempX = colshape(demeanedX[.,k],T)
		tempHx = rm_matvecmult(tempH, tempX)  // gives me a row vector indexed by t, with elements sum_s (H_ij,st * x_ijs)
		
		tempHx = tempHx # J(T,1,1)
		xHS_ij =  colshape(tempHx :* e,T*T)   // "e" is just the score. thus, this is the outer product of Hx and S.
		
		xHS_i = (collapse_i * xHS_ij)

		xHS_j = (collapse_j * xHS_ij)
		
		// compute xG objects that appear in the second terms in the bias.
		// tempX = colshape(demeanedX,T)#J(T,1,1)                                              // was: colshape(demeanedX[.,k],T)#J(T,1,1)
		
		//xG_ij = colshape(rm_matvecmult(tempG, tempX), T*T)	
		tempHx = colshape(H_ij :* demeanedX[.,k],T*T)
		xG_ij  = -tempHx - rm_transpose(tempHx)  // this does not give you Gx, but will give you equivalent sums when summed over i or over j.
		
		xG_i = (collapse_i  * xG_ij)
		xG_j = (collapse_j  * xG_ij)
		
		/*
		tempX = tempX#J(T,1,1)                                              // was: colshape(demeanedX[.,k],T)#J(T,1,1)
		xG_ij = colshape(rm_matvecmult(tempG, tempX), T*T)
		xG_i = (collapse_i  * xG_ij)
		xG_j = (collapse_j  * xG_ij)
		xG_ij
		*/
		
		xG_i = (collapse_i  * xG_ij)
		xG_j = (collapse_j  * xG_ij)

		// compute B1 and D1
		B1[.,k] = rm_matmult(xHS_i, H_i_inv) * fast_trace

		D1[.,k] = rm_matmult(xHS_j, H_j_inv) * fast_trace

		// compute B2 and D2		
		temp = rm_matmult(H_i_inv,xG_i)
		temp = rm_matmult(temp,H_i_inv)
		/*
		if (V_name != "") {
			temp_h = rm_matmult(temp,SSh_i) 		
			B2_h   = temp_h * fast_trace 
		}
		*/

		temp = rm_matmult(temp,SS_i) 
		B2[.,k]  = temp * fast_trace 

	    temp = rm_matmult(H_j_inv,xG_j)
	    temp = rm_matmult(temp,H_j_inv)
		/*
		if (V_name != "") {
			temp_h = rm_matmult(temp,SSh_j) 
			D2_h   = temp_h * fast_trace 
		}
		*/
		temp = rm_matmult(temp,SS_j) 
		D2[.,k]  = temp * fast_trace 
		
		B[k] = - (1/N_i) * sum( B1[.,k] ) + 1/(2*N_i) * sum( B2[.,k] )
		D[k] = - (1/N_j) * sum( D1[.,k] ) + 1/(2*N_j) * sum( D2[.,k] ) 
	}
		
	/*	
	"original bias correction"
	invsym(W) * (B / N_j_1 + D / N_i_1) 
	
	st_matrix(bias_name, invsym(W) * (B / N_j_1 + D / N_i_1) )	
	*/
	
	"bias corrections (to be subtracted from original coefficients)"
	invsym(W) * ((N_i / (N_i-1))*B / N_j_1 + (N_j / (N_j-1))*D / N_i_1 ) 		

	st_matrix(bias_name, invsym(W) * ((N_i / (N_i-1))*B / N_j_1 + (N_j / (N_j-1))*D / N_i_1) ) 
	
	/*
	if (V_name != "") {
		B_h  = - (1/N_i) * sum( B1 ) + 1/(2*N_i) * sum( B2_h )
		
		B_h2 = - (1/N_i) * sum( B1 ) + 1/(2*(N_i-1)) * sum( B2_h )
		
		D_h  = - (1/N_j) * sum( D1 ) + 1/(2*N_j) * sum( D2_h  ) 
		
		D_h2  = - (1/N_j) * sum( D1 ) + 1/(2*(N_j-1)) * sum( D2_h  ) 
		
		"bias_h"
		W^-1 * ((N_i / (N_i-1))*B_h / N_j_1 + (N_j / (N_j-1))*D_h / N_i_1) 
				
		st_matrix(bias_name+"_h", W^-1 * ((N_i / (N_i-1))*B_h / N_j_1 + (N_j / (N_j-1))*D_h / N_i_1) )	
	}
	*/
	
	st_matrix(W_name, W)
	
	st_matrix(B1_name, -sum(B1))	
	st_matrix(B2_name,  sum(B2)/2)
	st_matrix(B_name, B)
	
	st_matrix(D1_name, -sum(D1))
	st_matrix(D2_name,  sum(D2)/2)
	st_matrix(D_name, D)
}
end

mata: 
void twowayse_correction(string scalar demeaned_x, string scalar i_var, 
			             string scalar j_var, string scalar t_var, string scalar pair_var, string scalar y_var, string scalar lam, 
			             string scalar V_name, string scalar W_name, string scalar orig_sterr, 
			             | string scalar todo_var)
{

	vars = st_data(., (demeaned_x), todo_var)
	vars = st_data(., (i_var, j_var, t_var, pair_var, y_var, lam), todo_var), vars
	K   = cols(vars)-6
	
	T         = max(vars[.,3])
	NN_panels = max(vars[.,4]) 
	NNT = NN_panels * T
	
	ij = vars[.,1] + vars[.,2]/(2*max(vars[.,2]))
	ij = uniqrows(ij)#J(T,1,1)
	i  = floor(ij)
	j  = floor((ij-i):*(2*max(vars[.,2])) :+ 1e-6)
	
	t  = J(NN_panels,1,1)#(1..T)'
	it = i + t/(2*T)
	jt = j + t/(2*T)

	uniq_i = uniqrows(i)
	uniq_j = uniqrows(j)
	uniq_t = (1..T)'
	uniq_it = uniqrows(it)
	uniq_jt = uniqrows(jt)
	
	fast_trace = rowshape(I(T), T*T)
	
	N_i = rows(uniq_i)
	N_j = rows(uniq_j)
	
	N_j_1 = NN_panels/N_i
	N_i_1 = NN_panels/N_j
	
	// fill out panels with missing years
	index =  T :* (vars[.,4]:-1) :+ vars[.,3]
	tempvars = J(NNT, 3+K, 0)
	
	// fill out panels with missing years
	index =  T :* (vars[.,4]:-1) :+ vars[.,3]
	tempvars = J(NNT, 2+K, 0)
	tempvars[index,.] = editmissing(vars[., 5..(6+K)],0)
	
	y             = tempvars[.,1]
	lambda        = tempvars[.,2]
	demeanedX     = tempvars[.,3..(2+K)]
	
	e = y - lambda	
	
	// will use clustered standard errors here, which reduces to robust standard errors when T=1.
	W = (1/NN_panels) * ((lambda:*demeanedX)' * demeanedX) 
	
	// SS_ij is the OUTER product of scores
	SS_ij = (e # J(1,T,1)) :*(colshape(e,T)#J(T,1,1))
	SS_ij = colshape(SS_ij, T*T)
	
	// construct OMEGA^U
	// - demeanedX  is NNT x K (think of as NN T x K matrices, stacked vertically)
	// - demeanedX' is K x NNT (think of a NN K x T matrices, stacked horizontally)
	// - SSh_ij is NN x TT (think of as NN TT row vectors, stacked vertically)
	
	if(N_i*T+N_j*T<1000) {
		// sort of a fast xi 
		d_it = it :== (uniq_it')#J(NNT,1,1)
		d_jt = jt :== (uniq_jt')#J(NNT,1,1)
		d_ijt = (d_it, d_jt)
	
		// This can be shown to be equal to inv (sum_ij (d_ij' Lambda_ij d_ij))
		V_FE = invsym((d_ijt:*lambda)' * (d_ijt))

		// Next, we need d_ij V_FE d_ij' 
		//   - d_ijt  is NNT x F  (NN vertically stacked T x F matrices)
		//   - V_FE   is F x F    (NN stacked T x F matrices)
		//   - d_ijt' is F x NNT  (NN horizontally stacked F x T matrices)

		// d_ijt x V_F gives me another NNT x F matrix. 
		// - the elements are sum_f {d_ij1,f V_f1} sum_f {d_ij1,f V_f2} sum_f {d_ij1,f V_f3} ...
		//                    sum_f {d_ij2,f V_f1} sum_f {d_ij2,f V_f2} sum_f {d_ij2,f V_f3} ... 
		//                    sum_f {d_ij3,f V_f1} sum_f {d_ij3,f V_f2} sum_f {d_ij3,f V_f3} ... 
	
		// Now consider (d_ijt x V_F) :* d_ijt
		// - the elements would be sum_f {d_ij1,f V_f1} d_ij1,1 sum_f {d_ij1,f V_f2}d_ij1,2 sum_f {d_ij1,f V_f3} d_ij1,3 ...
		//                         sum_f {d_ij2,f V_f1} d_ij2,1 sum_f {d_ij2,f V_f2)d_ij2,2 sum_f {d_ij2,f V_f3} d_ij2,3 ... 
		//                         sum_f {d_ij3,f V_f1} d_ij3,1 sum_f {d_ij3,f V_f2)d_ij3,2 sum_f {d_ij3,f V_f3} d_ij3,3 ... 
	
		// Finding the rowsum then gives you the main diagonal.
		// To get the full outer product, pull apart the lefthand side so that it repeats.
		// For the right-hand side, reshape, pull apart, then reshape again.
		
		dVd  = J(NN_panels, T*T,0)
		dV    = d_ijt * V_FE
		dVd[.,(0..(T-1))*T :+ (1..(T))] = colshape(rowsum(dV :* d_ijt),T)
		
		/* // unnecessary; all off-diagonal terms in each of these TxT blocks are zero.
		if (T>=2) {
			for (u=2; u <= T; u++) {
				dVd[.,(0..(T-u))*T :+ (u..(T))] = colshape(rowsum(dV[selectindex(t:<(T-u+2)),.] :* d_ijt[selectindex(t:>(u-1)),.]),T-(u-1))	
	
				if(u==2) {
					idx1 = ((u-1)..(T-1))*T:+(u-1)
					idx2 = (u-2):+(u..T)
				}
				else{
					idx1 = idx1,((u-1)..(T-1))*T:+(u-1)
				idx2 = idx2,(u-2)*T:+(u..T)
				}
			}
			dVd[.,idx1] = dVd[.,idx2]  // this is the NN TxT blocks that lie aloneg the diagonal of d_ij V_FE d_ij'; each TxT block is laid out as a (TxT)x1 vector.
		}
		*/
	}
	else {
		alt_order     = order((i,j,t), (3,1,2))
				
		alt_inv_order = invorder(alt_order)     
		alt_i = i[alt_order]
		alt_j = j[alt_order]
		alt_lambda  = lambda[alt_order,.] 
				
		for (k=1; k <= T; k++) {
			t_index = (NN_panels*(k-1)+1..NN_panels*k)

			alt_uniq_i = uniqrows(alt_i[t_index])
			alt_uniq_j = uniqrows(alt_j[t_index])

			alt_d_ijt = alt_i[t_index,1] :== alt_uniq_i'#(J(NN_panels,1,1)), 
						alt_j[t_index,1] :== alt_uniq_j'#(J(NN_panels,1,1)) 
				
			V_FE_t   = invsym((alt_d_ijt:*alt_lambda[t_index])' * (alt_d_ijt))
		
			if (k == 1) {
				dVd_alt  = rowsum( (alt_d_ijt * V_FE_t ) :* alt_d_ijt)
			}
			else {
				dVd_alt  = dVd_alt \ rowsum( (alt_d_ijt * V_FE_t ) :* alt_d_ijt)
			}
		}
		dVd = dVd_alt[alt_inv_order,.]
		dVd = colshape(dVd :* (J(NN_panels,1,1) # I(T)), T*T)
	}
	
	// This gets us x_ij V_x x_ij'
	V_X = (1/NN_panels) * invsym(W)
	
	tempX = colshape(demeanedX, T*K)#J(T,1,1)
	tempX = colshape(tempX, K) 
	
	xVx = rowsum( ((demeanedX * V_X) # J(T,1,1)) :* tempX)
	xVx = colshape(xVx, T*T)
	
	Lambda_ij = lambda:*(J(NN_panels,1,1)#I(T))
	lev_correction = rm_matmult(colshape(Lambda_ij,T*T), xVx+dVd)
	lev_correction = rm_newtinv(colshape(fast_trace,T*T):-lev_correction,30,1e-12)
	
	SSh_ij = rm_matmult(lev_correction, SS_ij)

	tempX     = colshape (demeanedX, K*T)	
	X_reshape = J(1,K,1)#(1..T) :* K :+ (-(K-1)..0)#J(1,T,1)
	tempX = colshape(tempX[.,X_reshape],T)           // reshapes demeanedX so that it is now NNK x T - or NN K x T matrices vertically stacked.	
	
	tempO = rm_matvecmult(SSh_ij#J(K,1,1), tempX)    // gives me NN K x T matrices, stacked vertically
		
	// next I want to construct NN KxK outer products 
	tempO = (tempX#J(K,1,1)) :* colshape(colshape(tempO,T*K)#J(K,1,1),T) 
		// (I)  NN blocks; within each block, there are KK rows of T columns, with each k repeated K times
		// (II) NN blocks; within each block, there are KK rows of T columns, with 1..K repeated K times 
	
	tempO = rowsum(tempO)
	tempO = colsum(colshape(tempO, K*K))
		
	OMEGA = (1/(NN_panels)) * colshape(tempO,K)
	
	V     = (1/(NN_panels)) * invsym(W) * OMEGA * invsym(W)
	
	//  Stata uses a finite-sample correction to calculate SEs. See: 
	// 	- https://www.stata.com/meeting/13uk/nichols_crse.pdf.
	//  - https://www.stata.com/manuals/u20.pdf (p. 52)
	//	- https://www.stata.com/manuals/p_robust.pdf (p. 13)
	V = NN_panels / (NN_panels-1) * V
		
	"Adjusted SEs"
	sqrt(diagonal(V))
	st_matrix(V_name, V)

	//original VCV matrix
	tempO = rm_matvecmult(SS_ij#J(K,1,1), tempX)
	tempO = (tempX#J(K,1,1)) :* colshape(colshape(tempO,T*K)#J(K,1,1),T)
	tempO = rowsum(tempO)
	tempO = colsum(colshape(tempO, K*K))
	OMEGA = (1/(NN_panels)) * colshape(tempO,K)
	V_orig = (1/(NN_panels)) * invsym(W) * OMEGA * invsym(W)
	V_orig = NN_panels / (NN_panels-1) * V_orig

	//original SEs (for comparison)
	st_matrix(orig_sterr, (diagonal(sqrt(V_orig))))
}
end
