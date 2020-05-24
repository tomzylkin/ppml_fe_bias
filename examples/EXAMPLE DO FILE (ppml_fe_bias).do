*! Bias corrections for FE-PPML models with multi-way fixed effects
*!
*! Suggested citation: Weidner & Zylkin (2020): 
*! "Bias and Consistency in Three-way Gravity Models"
*! arXiv preprint arXiv:1909.01327

clear all
*global this_dir  "E:\Google drive\TOM LOCAL\Bias and Consistency FE PPML\Code\files to upload"
global this_dir   "C:\Users\ztom\Google Drive\TOM LOCAL\Bias and Consistency FE PPML\Code\files to upload"

cd "$this_dir"

cap set matsize 800
cap set matsize 11000
cap set maxvar 32000


// download data from github:
use "https://github.com/tomzylkin/ppml_fe_bias/blob/master/examples/PPMLFEBIAS_EXAMPLE_DATA.dta?raw=true" if category=="MANUF", clear // can change category to "NONMANUF" or "TOTAL"

// Manufacturing trade between 65 countries for the years 1986 - 2004, every four years.

// Sources:
//  - Trade: UN COMTRADE
// 	- FTAs:  NSF-Kellogg (Baier & Bergstrand) database
//  - "gravity" variables: CEPII (Head, Mayer, and Ries) gravity data


// estimate the average partial effect of an FTA on total trade
cap egen imp=group(isoimp)
cap egen exp=group(isoexp)
ppmlhdfe trade fta, a(imp#year exp#year imp#exp) cluster(imp#exp) d             // ssc install ppmlhdfe, replace

// save lambda values
predict lambda

// save coefficients
matrix beta = e(b)

// compute bias corrections
cap program drop ppml_fe_bias
ppml_fe_bias trade fta, i(exp) j(imp) t(year) lambda(lambda)

// create a results table
cap program drop ppml_fe_bias
ppml_fe_bias trade fta, i(exp) j(imp) t(year) lambda(lambda) beta(beta)

// twoway example
ppmlhdfe trade ln_distw contig colony comlang_off comleg fta, a(imp#year exp#year) cluster(imp#exp) d       

predict lambda_2way
matrix beta_2way = e(b)

cap program drop ppml_fe_bias
ppml_fe_bias trade ln_distw contig colony comlang_off comleg fta, i(exp) j(imp) t(year) lambda(lambda_2way) twoway beta(beta_2way)
