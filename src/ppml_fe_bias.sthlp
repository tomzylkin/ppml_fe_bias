{smcl}
{* *! version 1.0.1 24aug2016}{...}
{vieweralsosee "[R] xtpoisson" "help xtpoisson"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "ppml_panel_sg" "help ppml_panel_sg"}{...}
{vieweralsosee "ppmlhdfe" "help xtpqml"}{...}
{vieweralsosee "poi2hdfe" "help poi2hdfe"}{...}
{vieweralsosee "ppml" "help ppml"}{...}
{vieweralsosee "ge_gravity" "help ge_gravity"}{...}
{viewerjumpto "Syntax" "ppml_fe_bias##syntax"}{...}
{viewerjumpto "Description" "ppml_fe_bias##description"}{...}
{viewerjumpto "Main Options" "gppml_fe_bias##main_options"}{...}
{viewerjumpto "Background" "ppml_fe_bias##backgroup"}{...}
{viewerjumpto "Examples" "ppml_fe_bias##examples"}{...}
{viewerjumpto "Advisory" "ppml_fe_bias##advisory"}{...}
{viewerjumpto "Acknowledgements" "ppml_fe_bias##acknowledgements"}{...}
{viewerjumpto "References" "ppml_fe_bias##references"}{...}
{viewerjumpto "Citations" "ppml_fe_bias##citations"}{...}
{title:Title}

{p2colset 5 22 23 2}{...}
{p2col :{cmd:ppml_fe_bias} {hline 2}} Bias corrections for Poisson Pseudo-Maximum Likelihood (PPML) gravity models with two-way and three-way fixed effects.{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2}{cmd:ppml_fe_bias}
{depvar} [{indepvars}] 
{ifin}{cmd:,} {opt lambda(varname)} {opt i(exp_id)} {opt j(imp_id)} {opt t(time_id)} [{help ppml_fe_bias##options:options}] {p_end}

{p 8 8 2}{it: exp_id}, {it: imp_id}, and {it: time_id} are variables that respectively identify 
the origin, destination, and time period associated with each observation. "lambda" is an input for the conditional mean
of each observation. For more details, see the {browse "https://github.com/tomzylkin/ppml_fe_bias/blob/master/help%20file%20(ppml_fe_bias).pdf":online version of this help file}.{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:ppml_fe_bias}  implements analytical bias corrections described in Weidner & Zylkin (2020) for PPML "gravity"
regressions with two-way and three-way fixed effects, as are commonly used with international trade data and other types
of spatial flows. As shown by Weidner & Zylkin (2020), when the time dimension is fixed, the point estimates produced by
the three-way PPML gravity model have an asymptotic incidental parameter bias of order 1/N, where N is the number of countries,
and the cluster-robust variance estimator itself has a downward bias that is also of order 1/N. For the two-way PPML gravity model, only the standard errors are biased.
{p_end}

{marker main_options}{...}
{title:Main Options}

These options allow you to store results for the bias corrections, select the type of gravity model being used (two-way or three-way),
and stipulate whether an approximation should be used for the bias-corrected variance matrix.

{synoptset 20 tabbed}{...}
{synopt: {opt bias(name)}}Store bias corrections for coefficients as a Stata matrix.{p_end}

{synopt: {opt v(vame)}}Store bias-corrected variance matrix as a Stata matrix.{p_end}

{synopt: {opt w(name)}}Store the expected Hessian ("W") as a Stata matrix.{p_end}

{synopt: {opt b_term(name)}}Store the estimated "B"-component of the bias in the score for the three-way model that is associated with the origin-time fixed effects.{p_end}

{synopt: {opt d_term(name)}}Store the estimated "D"-component of the bias in the score for the three-way model that is associated with the destination-time fixed effects.{p_end}

{synopt: {opt beta(name)}}Pass the coefficients from a prior PPML estimation. When beta coefficients are provided,
a new results table will be produced complete with bias-corrected coefficients and standard
errors.{p_end}

{synopt: {opt twoway}}Pass the coefficients from a prior PPML estimation. When beta coefficients are provided,
a new results table will be produced complete with bias-corrected coefficients and standard
errors.

{synopt: {opt approx}}If “approx” is enabled, the bias correction for the variance will be computed using an approximation.
By default, this approximation is used whenever the number of origin-time and
destination-time fixed effects exceeds 1000 in order to facilitate computation and to avoid
running up against memory constraints.

{synopt: {opt exact}}Use an exact method for computing the bias-corrected variance, even if the number of origintime
and destination-time fixed effects exceeds 1000.

{synopt: {opt notable}}Suppress results table.

{marker background}{...}
{title:Background}

{pstd}Quick summary: {p_end}
{p 6 6 2} •  "Two-way" and "three-way" PPML gravity models are very popular for the estimation of international data and other types of spatial flows. 
The two-way model features origin-time and destination-time fixed effects. The three-way model then adds an origin-destination fixed effect.{break}
 •  Though three-way PPML gravity estimates are consistent, they suffer from an asymptotic bias due to the slow convergence of the origin-time and destination-time fixed effects. This bias causes confidence intervals to be incorrectly centered even in moderately large samples.{break}
 •  For a similar reason, cluster-robust standard errors that are typically used with three-way PPML are generally biased downward.{break}
 •  The two-way PPML gravity model gives estimates that are asymptotically unbiased but still suffers from downward-biased standard errors in much the same way as the three-way model.{break}
 •  This command implements analytical bias corrections that address each of these issues. These corrections are based on Taylor expansions described in Weidner & Zylkin (2020).{p_end}

{pstd}PPML gravity models with two-way and three-way fixed effects are very popular in the study of international trade
and other similar applications that involve bilateral flows (such as urban commuting or interregional migration).
The two-way PPML gravity model may be written as{p_end}

{p 8 8 2} y_ijt = exp[a_it + a_jt + x_ijt'b]w_ijt,{p_end}

{pstd}where y_ijt is a flow from origin i to destination j at time t, a_it and a_jt respectively are origin-time and destination-time fixed effects, b 
are the coefficients we want to estimate (typically variables that influence bilateral frictions, such as the distance between i and j), and 
w_ijt>=0 serves as an error term. The three-way model then adds a "country-pair" fixed effect a_ij:{p_end}

{p 8 8 2}y_ijt = exp[a_it + a_jt + a_ij + x_ijt'b]w_ijt.{p_end}

{pstd}In the latter model, the advantage of the third fixed effect a_ij is that it absorbs all time-invariant determinants of flows between i and j.
As discussed in Baier & Bergstrand (2007), the use of these fixed effects allows us to identify the elements of b based on time-variation in 
trade within pairs after first conditioning on a_it and a_jt. This approach is especially well suited for estimating the effects of bilateral
trade agreements and other similar policy variables and is currently recommended for this purpose by several leading references on gravity estimation
(e.g., Head & Mayer, 2014; Yotov, Piermartini, Monteiro, and Larch, 2016). {p_end}

{pstd}Weidner and Zylkin (2020)'s analysis identifies several econometric issues that arise with the three-way model, but their method for correcting the
standard errors of the three-way model also can be adapted to address a similar issue with the two-way model. Focusing on the three-way model for now,
it is important to recognize that the first-order conditions of PPML allow us to "profile out" the pair fixed effect a_ij from the model without biasing the scores of the other parameters
(see the {browse "https://github.com/tomzylkin/ppml_fe_bias/blob/master/help%20file%20(ppml_fe_bias).pdf":online version of this help file} for more details.)
Three-way PPML therefore has the nice property that it is consistent despite the large number of fixed effects, whereas other popular alternatives (OLS, Gamma PML) 
can be shown to be inconsistent in this setting. At the same time, the
two-way representation of the model also brings to mind the results from Fernández-Val & Weidner (2016) for the
asymptotic bias of two-way nonlinear models. The three-way PPML gravity model is a more complicated model
than the ones studied in Fernández-Val & Weidner (2016), but Weidner & Zylkin (2020) nonetheless demonstrate
that an analogous result occurs for three-way PPML whenever the time dimension is fixed.
Specifically, if I is the number of origins and J is the number of destinations,
then estimates for b in finite samples will have an asymptotic bias of the form{p_end}
	{p 8 8 2}(1/I) B* + (1/J) D*,{p_end}

{pstd}that is, a bias that vanishes only as both I and J -> infinity. Because the asymptotic standard error itself shrinks with 1/sqrt(IJ) as I and J become large, we have the discomfiting result
that our confidence intervals will be systematically off-center even in moderately large samples because of the large relative magnitude of the bias in relation to the standard error. 

{pstd}Weidner and Zylkin (2020) provide a more detailed derivation of the bias in b based on a second-order Taylor series expansion
of the expected profile score around around the correct values of the a_it and a_jt fixed effects. This expansion provides the basis
for the analytical bias correction performed by this command. In addition to this correction, {cmd:ppml_fe_bias} also addresses 
a related issue that affects the cluster-robust standard errors that are typically used with the three-way model. 
The latter problem is effectively a version of the more general result
that "heteroscedasticity-robust" standard errors tend to be downward-biased in small samples (cf., MacKinnon &
White, 1985; Imbens & Kolesar, 2016), only this problem is exacerbated here by the slow convergence of the fixed
effects, which causes the bias in the standard error to vanish at a much slower rate. A similar issue also arises in the 
two-way model (Egger & Staub, 2015; Jochmans,
2016; Pfaffermayr, 2019); thus, ppml_fe_bias has been programmed to provide bias-corrected standard errors
for two-way models as well as three-way models (by making use of the "twoway" option in the former case). The
b-coefficients obtained using the two-way PPML gravity model do not suffer from any asymptotic bias, however,
as originally shown by Fernández-Val & Weidner (2016).{p_end}

{marker examples}{...}
{title:Examples}

{pstd}These examples follow {browse "https://github.com/tomzylkin/ppml_fe_bias/blob/master/examples/EXAMPLE%20DO%20FILE%20(ppml_fe_bias).do":sample .do file} included along with this command. 
The data set used in this .do file
consists of a panel of 65 countries trading with one another over the years 1988-2004, using every 4 years. The
trade data uses aggregated trade flows from UN COMTRADE, with information on FTAs taken from the {browse "https://sites.nd.edu/jeffrey-bergstrand/database-on-economic-integration-agreements/": NSF-Kellogg database}
maintained by Scott Baier and Jeff Bergstrand and other covariates taken from the CEPII gravity
data set created by Head, Mayer, & Ries (2010). The computation of the PPML regression relies on the {browse "https://ideas.repec.org/c/boc/bocode/s458622.html":ppmlhdfe}
Stata command created by Correia, Guimarães, & Zylkin (2020).{p_end}

{pstd}{bf:Three-way example}. Using {cmd:ppmlhdfe}, the appropriate syntax for specifying a three-way gravity model
with exporter-time, importer-time, and exporter-importer fixed effects and standard errors that are clustered by pair is  {p_end}

{p 8 15 2}{cmd:ppmlhdfe trade fta, a(imp#year exp#year imp#exp) cluster(imp#exp) d}

{pstd}The "d" option is needed to facilitate obtaining values for the conditional mean of the dependent variable, which
{cmd:ppml_fe_bias} will need in order to construct the necessary expressions for the bias corrections. We next create
a new variable "lambda" containing the conditional mean from the regression as well as a matrix "beta" for the
estimated coefficient on fta. We then pass these results along with the data to {cmd:ppml_fe_bias}:{p_end}

{p 8 8 2}{cmd:predict lambda}{break}
{cmd:matrix beta = e(b)}{break}
{cmd:ppml_fe_bias trade fta, i(expcode) j(impcode) t(year) lambda(lambda) beta(beta)}{p_end}

{pstd}The resulting output shows there is an upward bias of about 0.006 in the estimated PPML coefficient for fta.
While the bias-corrected estimate is not overly different from the original uncorrected one (0.170 vs 0.178), it is
important to keep in mind this bias is about 22% of the estimated standard error, which is easily large enough to
make a meaningful difference for hypothesis testing in general cases. Consistent with the results of Weidner &
Zylkin (2020), the estimated standard error is also found to be biased: the bias-corrected standard error is about
20% larger than the uncorrected standard error, and the lower bound of the estimated 95% confidence interval
after both of these corrections are applied is substantially lower than what would be found otherwise (0.086 versus
0.109). While these results were obtained for a moderate number of countries (I=J=65), it is important to note that Weidner & Zylkin (2020)’s results imply that the magnitudes of
these biases are likely to depend significantly on the distribution of the data, even for ostensibly large data sets. Thus, it 
is recommended researchers implement these checks whenever feasible as a matter of good practice.{p_end}

{pstd}{bf:Two-way example}. Using the same example data set and .do file a typical two-way gravity regression would be

{p 8 8 2}{cmd:ppmlhdfe trade ln_distw contig colony comlang_off comleg fta, a(imp#year exp#year) cluster(imp#exp) d}{p_end}

{pstd}where we now include some covariates that would otherwise be absorbed by the exporter-importer fixed effect
from the three-way model, such as the log of bilateral distance and the sharing of a colonial relationship. The code
to compute bias-corrected standard errors and to present the results in a nicely formatted table would be{p_end}

{p 8 8 2}{cmd:predict lambda_2way}{break}
{cmd:matrix beta_2way = e(b)}{break}
{cmd:ppml_fe_bias trade ln_distw contig colony comlang_off comleg fta, i(expcode) j(impcode) t(year) lambda(lambda_2way) beta(beta_2way) twoway}{p_end}

{pstd}As the results show, the implied downward biases in the standard error for the two-way model are similar to the
corresponding bias found for fta coefficient in the three-way model, generally ranging between 20% and 26%.
The coefficients themselves are unbiased in this case, as shown by Fernández-Val & Weidner (2016).{p_end}

{pstd}{bf:Storing results}.So long as an input is given for the “beta” matrix, {cmd:ppml_fe_bias} will store results for the bias-corrected
coefficients and variance matrix using {cmd:ereturn post}. This allows users to access these results using
the post-estimation operators {cmd:e(b)} and {cmd:e(V)}. It also makes it possible to use {cmd:ppml_fe_bias} in conjunction with
standard table-formatting commands such as {cmd:estout} or {cmd:estimates table}. {p_end}

{pstd}{bf:Approximation for the adjusted variance}. Depending on the size of the data, it may be necessary to use an
approximation method to compute the necessary bias correction for the cluster-robust variance matrix. Because
the details behind this approximation are mathematically complex, they have been left for the 
{browse "https://github.com/tomzylkin/ppml_fe_bias/blob/master/help%20file%20(ppml_fe_bias).pdf":online version of this help file}. In brief, the method
involves constructing a convergent sequence for the variance correction based on alternating projections and taking the result after the first few iterations.
Doing so manages to avoid any large matrix operations that would otherwise be necessary and generally yields a reasonable approximation that still leads to significantly
improved inferences. For researchers who would prefer not to use this approximation, the "exact" option may be used to force {cmd:ppml_fe_bias} to calculate the adjusted variance directly.{p_end}

{marker other}{...}
{title:Advisory}

{pstd}This is an advanced technique that requires a basic understanding of PPML estimation and of the three-way
gravity model. I would recommend either Yotov, Piermartini, Monteiro, & Larch (2016) or Larch, Wanner, Yotov,
& Zylkin (2019) for further reading on these topics. For essential reading on two-way gravity estimation, see
Head & Mayer (2014).{p_end}

{pstd}This is version 1.0 of this command. Depending
on interest, future versions of this command could add further options for multi-way clustering, multi-industry
models, and/or dynamic models. If you believe you have found an error that can be
replicated, or have other suggestions for improvements, please feel free to {browse "mailto:tomzylkin@gmail.com":contact me}.{p_end}

{marker contact}{...}
{title:Author}

{pstd}Thomas Zylkin{break}
Department of Economics, Robins School of Business{break}
University of Richmond{break}
Email: {browse "mailto:tomzylkin@gmail.com":tomzylkin@gmail.com}
{p_end}

{marker citation}{...}
{title:Suggested Citation}

If you are using this command in your research I would appreciate if you would cite

{pstd}• Weidner, Martin and Thomas Zylkin. “Bias and Consistency in Three-way Gravity Models." arXiv preprint
arXiv:1909.01327 (2020).{p_end}

The code used in this command implements the bias corrections described in Section 3 of our paper for the threeway
model. The bias correction for the variance of the two-way model is discussed in the appendix.

{marker acknowledgements}{...}
{title:Acknowledgements}

{pstd} The computations needed to construct these bias corrections have been greatly facilitated by three outstanding
Stata packages: {browse "https://ideas.repec.org/c/boc/bocode/s457888.html":rowmat_utils} by Matthew J. Baker, {browse "https://ideas.repec.org/c/boc/bocode/s458514.html": gtools} 
by Mauricio Caceres Bravo, and {browse "https://ideas.repec.org/c/boc/bocode/s457985.html":hdfe} by Sergio Correia.
The code for displaying results uses the {browse "https://econpapers.repec.org/scripts/search.pf?ft=frmttable":fmttable} command by John Luke Gallup. Thanks to each of these
creators for these incredibly useful tools!.{p_end}

{marker further_reading}{...}
{title:Further Reading} 

{pstd}• Gravity estimation: Head and Mayer (2014); Yotov, Piermartini, Monteiro, & Larch (2016){p_end}

{pstd}• Bias corrections for two-way models: Fernández-Val & Weidner (2016){p_end}

{pstd}• Bias corrections for other three-way models aside from PPML: Fernández-Val & Weidner (2018); Hinz, Stammann, and Wanner (2019) {p_end}

{pstd}• Correia, Guimarães, and Zylkin (2020); Larch, Wanner, Yotov, and Zylkin (2019); Stammann (2018) {p_end}


{marker references}{...}
{title:References}

{phang}
Baier, S. L. & Bergstrand, J. H. (2007), “Do Free Trade Agreements actually Increase Members’ International
Trade?”, Journal of International Economics 71(1), 72–95.{p_end}

{phang}
Correia, S., Guimarães, P., & Zylkin, T. (2020), “Fast Poisson Estimation with High-dimensional
Fixed Effects”, Stata Journal 20(1), 95-115.{p_end}

{phang}
Egger, P. H. & Staub, K. E. (2015), “GLM estimation of trade gravity models with fixed effects”, Empirical
Economics 50(1), 137–175.{p_end}

{phang}
Fernández-Val, I. & Weidner, M. (2016), “Individual and Time Effects in Nonlinear Panel Models with Large N,
T”, Journal of Econometrics 192(1), 291–312.{p_end}

{phang}
Fernández-Val, I. &Weidner, M. (2018), “Fixed Effect Estimation of Large T Panel Data Models”, Annual Review
of Economics.{p_end}

{phang}
Head, K. & Mayer, T. (2014), “Gravity Equations: Workhorse, Toolkit, and Cookbook”, Handbook of International
Economics 4, 131–196.{p_end}

{phang}
Head, K., Mayer, T., & Ries, J. (2010), “The erosion of colonial trade linkages after independence”, Journal of
International Economics 81(1), 1–14.{p_end}

{phang}
Hinz, J., Stammann, A., & Wanner, J. (2019), “Persistent Zeros: The Extensive Margin of Trade”, Unpublished
manuscript.{p_end}

{phang}
Imbens, G. W. & Kolesar, M. (2016), “Robust standard errors in small samples: Some practical advice”, Review
of Economics and Statistics 98(4), 701–712.{p_end}

{phang}
Jochmans, K. (2017), “Two-Way Models for Gravity”, Review of Economics and Statistics 99(3), 478-485.{p_end}

{phang}
Larch, M., Wanner, J., Yotov, Y. V., & Zylkin, T. (2019), “Currency Unions and Trade: A PPML Re-assessment
with High-dimensional Fixed Effects”, Oxford Bulletin of Economics and Statistics 81(3), 487–510.{p_end}

{phang}
MacKinnon, J. G. & White, H. (1985), “Some Heteroskedasticity-consistent Covariance Matrix Estimators with
Improved Finite Sample Properties”, Journal of Econometrics 29(3), 305–325.{p_end}

{phang}
Pfaffermayr, M. (2019), “Gravity models, PPML estimation and the bias of the robust standard errors”, Applied
Economics Letters 26(18), pp. 1–5.{p_end}

{phang}
Stammann, A. (2018), “Fast and feasible estimation of generalized linear models with high-dimensional k-way
fixed effects”, arXiv preprint arXiv:1707.01815.{p_end}

{phang}
Weidner, M. & Zylkin, T. (2020), “Bias and Consistency in Three-Way Gravity Models”, arXiv preprint
arXiv:1909.01327.{p_end}

{phang}
Yotov, Y. V., Piermartini, R., Monteiro, J.-A., & Larch, M. (2016), “An Advanced Guide to Trade Policy Analysis:
The Structural Gravity Model”, World Trade Organization, Geneva.{p_end}

