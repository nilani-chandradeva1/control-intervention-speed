
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Modelling the influence of the duration of effectiveness and time taken to implement vector control interventions for vector-borne diseases. 

Authors: Nilani Chandradeva*, Andrew C. Glover, Charles Whittaker, Joseph D. Challenger, Hannah C. Slater**, Ellie Sherrard-Smith** & Thomas S. Churcher**

Corresponding author*: Nilani Chandradeva (nilani.chandradeva@gmail.com)
**: Joint senior authors

The scripts in this repository are for running the analysis presented in the titled manuscript. 

The scripts for this analysis are found in the following folders : `1.simple-seirs-sei` and `2.ivm-stagger-distr`. Odin models within the `inst` folder are called from scripts in these two folders.
Within each folder, scripts are numbered in the order they should be run in. Output generated from each analysis are saved in the `output` folder and will be called from this folder if required later on. All figures, made using scripts in the `figures` folder, are saved in `plots`. This work builds on previously conducted work using the ivRmectin package (https://github.com/mrc-ide/ivrmectin), described below.

<!-- badges: start -->

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
<!-- badges: end -->

*ivRmectin* enables users to simulate models of malaria transmission
with an explicit focus on exploring the impact of different
endectocide-based interventions such as ivermectin. This builds on
previous work by Joel Hellewell, OJ Watson, Hannah Slater, Juliette
Unwin, Rich Fitzjohn and Hannah Slater who developed a deterministic
version of the Imperial College Malaria Model (for more information
about the original model which is stochastic and individual based, see
*Griffin et al., 2010 PLOS Medicine*:
<https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1000324>).

Subsequent work led by Hannah Slater and the rest of the team extended
this model to model the impact of endectocides such as ivermectin on
malaria epidemiology in endemic settings. This model has been previously
published and more information can be found in *Slater et al., 2014,
Journal of Infectious Diseases*
(<https://academic.oup.com/jid/article/210/12/1972/2908666>) and *Slater
et al., 2020, Lancet Infectious Diseases*
(<https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(19)30633-4/fulltext>).
It is this model that *ivRmectin* contains, and allows users to run.
