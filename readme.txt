This directory contains the necessary MATLAB code and files to replicate the simulation results in 
"Large Skew-t Copula Models and Asymmetric Dependence in Intraday Equity Returns"
Authors: Lin Deng, Michael Smith, Worapree Maneesoonthorn

MATLAB files:
'ACCop_VB_example.m'	Replicates the simulation results in Section 4 of the paper
'ACCop_VB_example_2.m'	Code for estimating currency exchange rate dependence using AC skew-t copula. 
			Note that this example is provided to illustrate how the code may be used in a real data setting, but is not covered in the paper.
'Generate_CopulaData.m'	Code for generating simulated data from the AC skew-t copula density

Directories:
'Data'		Folder containing simulated data for d=5,30 (used in the simulation exercise in the paper) and the exchange rate data used in 'ACCop_VB_example_2.m'. 
'Distribution'	Folder containing functions related to the evaluation of the AC skew-t distribution and copula
'misc'		Folder containing miscellaneous functions for saving and plotting outputs
'Results'	Folder containing simulation results, as presented in Section 4 of the paper
'VB_fun'	Folder containing functions related to the Variational Bayes estimation of the AC skew-t copula.
		Key functions in 'VB_fun': 	
			'vb_st_copula_opt_b.m' 	  estimates \lambda, the variational parameters of the skew-t copula and its special cases (the AC skew-normal, Gaussian and t copula)
			'summary_stc_vb.m' 	  calculates the summaries of the posterior approximation of the AC model parameters
			'sim_posterior_stc_vb.m'  simulates the copula data from the predictive posterior distribution
			'gradLogPost_TraceGrad01' computes the analytical gradient of the log-posterior given in Table 5 of the manuscript.

