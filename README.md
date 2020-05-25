# `ppml_fe_bias`: bias corrections for two-way and three-way fixed effects PPML models

- Current version: `1.0 24feb2020`
- Jump to: [`citation`](#citation) [`example output`](#example-output) [`install`](#installation)
- Also see: [Help file](https://github.com/tomzylkin/ppml_fe_bias/blob/master/help%20file%20(ppml_fe_bias).pdf) | [Examples](https://github.com/tomzylkin/ppml_fe_bias/blob/master/examples/EXAMPLE%20DO%20FILE%20(ppml_fe_bias).do) | [`ado file`](https://github.com/tomzylkin/ppml_fe_bias/blob/master/src/ppml_fe_bias.ado) | [`ppmlhdfe`](https://github.com/sergiocorreia/ppmlhdfe) 

**ppml_fe_bias** is a Stata package that implements analytical bias corrections described in [Weidner and Zylkin (2020)](https://arxiv.org/pdf/1909.01327.pdf) for PPML “gravity” regressions with two-way and three-way fixed effects, as are commonly used with international trade data and other types of spatial flows. As shown in [Weidner and Zylkin (2020)](https://arxiv.org/pdf/1909.01327.pdf), when the time dimension is fixed, the point estimates produced by the three-way PPML gravity model have an asymptotic incidental parameter bias of order `1/N`, where `N` is the number of countries, and the cluster-robust sandwich estimator that is typically used for inference itself has a downward bias that is also of order `1/N`. For the two-way PPML gravity model, only the standard errors are biased.

## Example output

Three-way fixed effects example

<p align="center"><img src="https://github.com/tomzylkin/ppml_fe_bias/blob/master/examples/figures/example%20output%20(3%20way).png?raw=true" alt="example output"/></p>

Two-way fixed effects example

<p align="center"><img src="https://github.com/tomzylkin/ppml_fe_bias/blob/master/examples/figures/example%20output%20(2way).png?raw=true" alt="example output (2 way)"/></p>

These examples follow the included [.do file and data set](https://github.com/tomzylkin/ppml_fe_bias/blob/master/examples/EXAMPLE%20DO%20FILE%20(ppml_fe_bias).do). Note that only the standard errors are biased in the case of the two-way model. In the three-way model, the point estimates also suffer from some bias. 

## Citation

Weidner, Martin & Thomas Zylkin: “Bias and Consistency in Three-way Gravity Models”, 2020; <a href='https://arxiv.org/pdf/1909.01327'>arXiv:1909.01327</a>.

## Installation

`ppml_fe_bias` requires the latest versions of [`hdfe`](https://ideas.repec.org/c/boc/bocode/s457985.html), [`gtools`](https://gtools.readthedocs.io/en/latest/), [`rowmat_utils`](https://ideas.repec.org/c/boc/bocode/s457888.html), and [`frmttable`](https://ideas.repec.org/c/boc/bocode/s375201.html) 

An SSC version that can be downloaded from the SSC repository within Stata will be available soon.

For now, you can install the latest version from Github:

```stata
cap ado uninstall ppml_fe_bias
net install ppml_fe_bias, from("https://raw.githubusercontent.com/tomzylkin/ppml_fe_bias/master/src") replace
```

After installing, you can then access a help file with more information in the standard way by typing

```stata
help ppml_fe_bias 
```

An online version of this help file with additional details is also available [here](https://github.com/tomzylkin/ppml_fe_bias/blob/master/help%20file%20(ppml_fe_bias).pdf). 


