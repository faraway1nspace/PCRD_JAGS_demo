# PCRD_JAGS_demo
A full-capture Hierarchical Bayesian model of Pollock's Closed Robust Design and application to dolphins

This contains a tutorial of several Robust Design (RD) Bayesian Mark-Recapture models, as presented in the paper <i>"Rankin RW, Nicholson K, Allen SJ, Krützen M, Bejder L, Pollock KH. 2016. A full-capture Hierarchical Bayesian model of Pollock's Closed Robust Design and application to dolphins. Frontiers in Marine Science, in Press".</i> This tutorial also uses (annonamized) data directly from the parent publication: <i>"Nicholson K, Bejder L, Allen SJ, Krützen M, and Pollock KH. 2012. Abundance, survival and temporary emigration of bottlenose dolphins (Tursiops sp.) off Useless Loop in the western gulf of Shark Bay, Western Australia. Marine and Freshwater Research 63:1059–1068."</i> Please cite both these studies if using this code and data (See BIBTEX citations at the end of this read me). 

The RD models are useful for the analysis of marked animals, estimating their abundance and demographic parameters, such as survival and movement. Importantly, the Hierarchical Bayesian (HB) model creates a flexible modelling environment to address several perennial issues encountered in Frequentist models 'fixed-effects' models. For example, HB results can be very similar to 'model-averaging' estimates. HB can also trivially model individual heterogeneity in detection probabilities (the bane of most mark-recapture studies).

<b> Dependencies </b>

This demo assumes the user has both R and JAGS installed, as well as the R package "rjags" to communicate between the two. See http://cran.r-project.org/ and http://mcmc-jags.sourceforge.net/. Most Linux distros have R/JAGS binaries available directly in standard repos. For example, Ubuntu users can simply type:
> sudo apt-get update

> sudo apt-get install r-base-core r-base-dev jags

Once in R, type: 

> install.package("rjags")

<b> Files in the tutorial </b>

There are 3 tutorials in 3 R files. Each tutorial presents a slightly different Bayesian model:
 - R_PCRD_JAGS_firstcapt_fixedeff.R 
 - R_PCRD_JAGS_fullcapt_fixedeff.R  
 - R_PCRD_JAGS_hierarchical.R (<b>USERS SHOULD START HERE!</b> see the animation about hyperpriors) 

The first file <i>R_PCRD_JAGS_firstcapt_fixedeff.R</i> is a tutorial for Pollock's Closed Robust Design, conditioning on first capture and is a "fixed-effects" version. The 2nd file <i>R_PCRD_JAGS_fullcapt_fixedeff.R</i> is also a fixed-effects version, but models the entire capture-history, thereby facilitating inference about recruitment processes (births). The 3rd file <i>R_PCRD_JAGS_hierarchical.R</i>is a Hierarchical Bayesian version, which will be of interest to most users: it can accommodate individual heterogeneity, full-capture modelling, and employs a peculiar type of "hyperprior" to induce a desirable shrinkage between fully time-varying parameters and time-constant parameters. In the companion paper, the authors argue that this hyperprior specification results in a type of "model-averaging" somewhat similar to the AICc-averaging that is ubiquitous in Mark-Recapture studies (thanks to the popularity of Program MARK). 

Other files are:
 - mark_capture_histories.inp
 - R_PCRD_JAGS_SOURCE.R
 - various ".JAG" files (model syntax) 

The .inp file is a MARK data file containing capture histories for 5 years of photo-ID studies of bottlenose dolphins in the western gulf of Shark Bay, as presented originally in "Nicholson et al. 2012." Please cite Nicholson when using this data. The SOURCE file contains some handy auxiliary functions to make it easier to initialize the JAGS models (especially generating sensible initial values for the latent-states of the HMM). The .JAG files are automatically generated when a user runs the above R tutorial scripts, and are not strictly necessary; they have the raw JAGS syntax for each JAGS model.

<b> Getting Started </b>

After cloning the appropriate files, users can re-analyze the data from Rankin et al. (2016) and Nicholson et al. (2012) by simply stepping through the R tutorials. Users can tweak priors and model assumptions by customizing the JAGS model syntax (see the R files).

<b> First-capture vs. Full-capture </b>

Users interested in either a) individual heterogeneity, or b) estimating recruitment processes should use the full-capture models. An advantage of the first-capture model is that users can model initial capture probabilities as being different from subsequent recapture probabilities (i.e., the p!=c specification in Program MARK). WARNING: users interested in the first-capture conditioning should familiarize themselves with the WinBUGS 'zeros trick'for exotic Likelihood functions (see http://users.aims.ac.za/~mackay/BUGS/Manuals/Tricks.html). Our First-capture models use a conditional multinomial distribution, the evaluation of which makes the JAGS Syntax somewhat more daunting for new users.

<b> Customizing For Other Studies: The Hierarchical Model Hyper-Priors </b>

Users interested in the Hierarchical Bayesian (HB) model will likely want to change the "hyperparameters" that govern the amount of shrinkage between time-varying and time-constant parametrizations. The hyperparameters are passed to jags as 'data', called "pr.tauphi","pr.taug1","pr.taug2","pr.taupdmu","pr.taupd2nd","pr.taueps" in the JAGS syntax.  Each of these hyperparameters govern a hyperprior distribution called the Scaled Half Student-T, with hyperparameters "mu","tau", and "nu". The half-T has favourable shrinkage-inducing properties, in that when there is little evidence in the data for widely-varying parameters, it shrinks a dispersion parameter to its prior expectation (near 0) and enforces a near time-constant parametrization for a particular parameter; conversely, lots of evidence in the data can result in large dispersion parameters and fully-time varying parameterizations, close to MLEs. Intermediate values for the Student T hyperparameters seem to result in parameter estimates similar to AICc model-averaged estimates.

See the animation. On the LEFT, there are different hyperpriors for sigma, each specified with different hyparameters. On the RIGHT, are the consequences of each hyperprior on individual p estimates. Under low-sample sizes, all the p's are "shrunk" to their overall group mean around 0.5. With more data, the p's approach their MLEs. This is parsimony: simple models under sparse data, and more complex models with larger sample sizes. This is often the goal of "model-averaging", but NOT for frequentist ML-based estimation. 

[[https://github.com/faraway1nspace/PCRD_JAGS_demo/blob/master/img/HalfTdemo.gif|alt=studentTdemo]]

For example in JAGS, the "pr.tauphi" hyperparameters will govern the disribution on "sigma.phi": if "sigma.phi" is close to zero, the survival "phi" will be nearly constant over time, and shrunk to its mean "phi.mu". Conversely, when "sigma.phi" is not zero, then "phi" will be fully time-varying, with individual phi values deviating strongly from their mean "phi.mu", and perhaps close to their individual Maximum Likelihood Estimates. The amount of shrinkage is controlled partially by the i) sample-size/amount of data, and ii) by the hyperparameters pr.tauphi=c(mean=0,tau=1/0.3^2,nu=2). Large values of 'tau' set the (approximate) prior expectation of 'sigma.phi', while 'nu' controls the tail of the distributions. Small values of 'nu' (2-3) allow the data more lee-way in deciding the value of 'sigma.phi'; large values of 'nu' make the hyper-prior more informative. A special case to avoid is nu=1, which is the half-Cauchy distribution, and allows too much dispersion for a random-variable that is between transformed between the logit and probability scales.

<b>We recommend tau >= 0.3^-2 and nu => 2 as good starting points for the half Student-t hyper-parameters (especially detection probabilities and gamma-double-prime)</b>. Such values will result in parameter values that are close to the MLEs (if there is a lot data), but will still be slighly shrunk towards their time-constant values under weak/sparse data. However, for those parameters that are more difficult to identify and estimate, even with a lot of data, like survival and gamma-prime, it may be desireable to enforce more shrinkage. For example, in the bottlenose dolphin example, we had a prior expectation that there would be little year-to-year variation in apparent survival (survival should only vary 1-2% per year), and so we effectively coerced a near time-constant phi with the hyperparameters tau = 1/0.2^2, nu=13. 

The above are hyperpriors, and therefore should be motivated by one's prior knowledge. Notice that in the hierarchical model, we place informative hyperpriors on the <i>dispersion</i> parameters, but admit little prior knowledge about the actual <i>means</i> of parameters.  

<b> Initializing Bayesian RD Capture-Recapture models </b>

The Capture-Recapture models are possible in JAGS because we can re-parameterize the Capture-Recapture process as a state-space model (or Hidden Markov Model). This means we have a stochastic time-series of 'latent states' z={<i>dead, inside, outside</i>} that are random variables just like any other random variable in the model. Unforunately, it means that we must provide initial values of z for JAGS, in order to start the MCMC chains. This is a pain in the posterior (!) for the novice user: ergo, we have provided a handy forwards-messaging/backwards-sampling algorithm to estimate initial values of these latent states z. See the "R_PCRD_JAGS_SOURCE.R".

<b>Citating: Bibtex</b>

If you use or modify these codes, please cite JAGS, R, and the companion paper below Rankin et al. 2016. If you use the bottlenose dolphin data from Useless Loop, please cite Nicholson et al. 2012. You can import the following into your reference manager:


@article{rankin_full-capture_2016,
	title = {A full-capture {Hierarchical} {Bayesian} model of {Pollock}'s {Closed} {Robust} {Design} and application to dolphins},
	journal = {Frontiers in Marine Science},
	author = {Rankin, Robert W. and Nicholson, Krista E. and Allen, Simon J. and Krützen, Michael and Bejder, Lars and Pollock, Kenneth H.},
	year = {2016}
}

@article{nicholson_abundance_2012,
	title = {Abundance, survival and temporary emigration of bottlenose dolphins (\textit{{Tursiops} sp.}) off {Useless} {Loop} in the western gulf of {Shark} {Bay}, {Western} {Australia}},
	volume = {63},
	doi = {10.1071/MF12210},
	number = {11},
	journal = {Marine and Freshwater Research},
	author = {Nicholson, Krista and Bejder, Lars and Allen, Simon J. and Krützen, Michael and Pollock, Kenneth H.},
	year = {2012},
	pages = {1059--1068}
}

This project is possible because of the Shark Bay Dolphin Innovation Project. See http://www.sharkbaydolphins.org/?page_id=90.
