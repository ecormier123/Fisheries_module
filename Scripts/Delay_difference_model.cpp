// SSModel with mods (SS2v4):
//   * clappers as obs eq (fit to Lt)
//   * random walk on mt, sd = sigma_m
//   * random walk on Rt (Ricker not useful), sd = sigma_phi
//   * obs eq for Ct following Thorson et al. (2013, CJFAS), sd = sigma_C
// Summary:
//   * 3 proc eq: biomass, recruits and nat mort (*natural mortality rate)
//   * 4 obs eq: clappers (* L: survey estimates of the number of scallops that died due to natural causes), survey commercial (* I: survey index of the biomass of commercial size scallops), survey recruits (* R: survey index of the biomass of recruitment size scallops), catch (* C: commercial catch reported at the end of the year)

#include <TMB.hpp>

template <class Type> 
Type sqr(Type x){
	return x*x;
}

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x,eps,Type(0.01)*pow(x-eps,2),Type(0));
  return CppAD::CondExpGe(x,eps,x,eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() () {
	
	//----------------------------------------------------------------------------
	// Data
	//----------------------------------------------------------------------------
    //(* Directly observable variable)
	DATA_VECTOR(I); // response: survey index commercial size (dim NY)
	DATA_VECTOR(IR); // response: survey index recruit size (dim NY)
	DATA_VECTOR(L); // response: survey estimate clappers (dim NY)
	DATA_VECTOR(C); // response: catch index (dim NY)
	DATA_VECTOR(g); // covariate: growth rate commercial size (dim NY)
	DATA_VECTOR(gR); // covariate: growth rate recruit size (dim NY)
	DATA_VECTOR(N); // covariate: survey estimate of live animals (dim NY)
	DATA_VECTOR(mobs); // covariate: survey estimate of ratio of dead to total number of live, estiamte of mortality.
	//Fit the alternative model without the fishing effort equation(9)
	//Taking in "enable_catcheq" as an indicator: 0 means excluding equation(9), 1 means including equation(9).
	DATA_INTEGER(enable_catcheq);
  // Do we want to incluce the popcorn model (1) or use the simipler alternative model (0).
	DATA_INTEGER(popcorn);
	//----------------------------------------------------------------------------
	// Parameters
	//----------------------------------------------------------------------------

	PARAMETER(log_sigma_tau); // proc sd biomass
	PARAMETER(log_sigma_phi); // proc sd recruits
	PARAMETER(log_sigma_m); // proc sd nat mort
	PARAMETER(log_sigma_epsilon); // obs sd survey indices commercial size
	PARAMETER(log_sigma_upsilon); // obs sd survey indices recruits
	PARAMETER(log_sigma_kappa); // obs sd clappers
	PARAMETER(log_sigma_C); // obs sd effort dynamic for catches
	PARAMETER(log_q_I); // catchability for commercial size
	PARAMETER(log_q_R); // catchability for recruits
	PARAMETER(log_S); // dissolution rate in clappers eq
	PARAMETER(log_a); // in effort dynamic catch eq
	PARAMETER(log_chi); // in effort dynamic catch eq

	PARAMETER_VECTOR(log_B); // randeff: biomass (dim NY)
	PARAMETER_VECTOR(log_R); // randeff: recruits (dim NY)
	PARAMETER_VECTOR(log_m); // randeff: nat mort (dim NY)
	

	//----------------------------------------------------------------------------
	// Setup
	//----------------------------------------------------------------------------

	int NY = I.size(); // number of years

	Type sigma_tau = exp(log_sigma_tau);
	Type sigma_phi = exp(log_sigma_phi);
	Type sigma_m = exp(log_sigma_m);
	Type sigma_epsilon = exp(log_sigma_epsilon);
	Type sigma_upsilon = exp(log_sigma_upsilon);
	Type sigma_kappa = exp(log_sigma_kappa);
	Type sigma_C = exp(log_sigma_C);

	Type q_I = exp(log_q_I);
	Type q_R = exp(log_q_R);
	//Type q_M = exp(log_q_M);
	Type S = exp(log_S);
	Type a = exp(log_a);
	Type chi = exp(log_chi);
	
	vector<Type> B = exp(log_B);
	vector<Type> R = exp(log_R);
	vector<Type> m = exp(log_m);
	vector<Type> nll_comp(7); // nll components, break down contributions
	nll_comp.setZero(); // initialize 

	//----------------------------------------------------------------------------
	// Proc eq
	//----------------------------------------------------------------------------

	// RW for recruits
	for (int t=1; t<NY; t++){ // R[0] free, yet predicted
	  	Type mean_proc_R = R[t-1];
	  	nll_comp[0] -= dnorm(log_R[t], log(mean_proc_R)-sqr(sigma_phi)/2.0, sigma_phi, true);
	  	// SIMULATE{
	  	//   log_R[t] = rnorm(log(mean_proc_R)-sqr(sigma_phi)/2.0, sigma_phi);
	  	// }
	}
	
	// RW for nat mort
	for (int t=1; t<NY; t++){ // m[0] free, yet predicted
	  	Type mean_proc_m = m[t-1];
	  	nll_comp[1] -= dnorm(log_m[t], log(mean_proc_m)-sqr(sigma_m)/2.0, sigma_m, true);
	  	// SIMULATE{
	  	//   log_m[t] = rnorm(log(mean_proc_m)-sqr(sigma_m)/2.0, sigma_m);
	  	// }
	}
	// biomass
	for (int t=1; t<NY; t++){ // B[0] free, yet predicted
	  Type mean_proc_B = exp(-m[t])*g[t-1]*(B[t-1]-C[t-1])+exp(-m[t])*gR[t-1]*R[t-1];
	  nll_comp[2] -= dnorm(log_B[t], log(mean_proc_B)-sqr(sigma_tau)/2.0, sigma_tau, true);
	}
	// Remove the penalty term option (from TMB 2.0 course, idea from Anders), leading to some weird results!
	// Idea to ensure function returns a positive result (i.e. C < B) for simulations.
	//Type pen=0;
	// for (int t=1; t<NY; t++){ // B[0] free, yet predicted
	//   	Type mean_proc_B = exp(-m[t])*g[t-1]*posfun((B[t-1]-C[t-1]),Type(1.0e-6),pen)+exp(-m[t])*gR[t-1]*R[t-1];
	//   	nll_comp[2] -= dnorm(log_B[t], log(mean_proc_B)-sqr(sigma_tau)/2.0, sigma_tau, true)+pen;
	//   	// SIMULATE{
	//   	//   log_B[t] = rnorm(log(mean_proc_B)-sqr(sigma_tau)/2.0, sigma_tau);
	//   	// }
	// }
	
	//----------------------------------------------------------------------------
	// Obs eq
	//----------------------------------------------------------------------------

	// survey indices for commercial size
	for (int t=0; t<NY; t++){
		Type mean_obs_I = q_I*B[t];
  	nll_comp[3] -= dnorm(log(I[t]), log(mean_obs_I)-sqr(sigma_epsilon)/2.0, sigma_epsilon, true);
  	// SIMULATE{
  	//   I[t] = exp(rnorm(log(mean_obs_I)-sqr(sigma_epsilon)/2.0, sigma_epsilon));
  	// }
	}

	// survey indices for recruits
	for (int t=0; t<NY; t++){
		Type mean_obs_R = q_R*R[t];
  	nll_comp[4] -= dnorm(log(IR[t]), log(mean_obs_R)-sqr(sigma_upsilon)/2.0, sigma_upsilon, true);
  	// SIMULATE{
  	//   IR[t] = exp(rnorm(log(mean_obs_R)-sqr(sigma_upsilon)/2.0, sigma_upsilon));
  	// }
	}

	// clappers, do we run the popcorn model or leave it out...
	if(popcorn == 1)
	{
	  Type mean_obs_L = m[0]*S*N[0]; // initial one different
	  nll_comp[5] -= dnorm(log(L[0]), log(mean_obs_L)-sqr(sigma_kappa)/2.0, sigma_kappa, true);
	  //SIMULATE{
	 //   L[0] = exp(rnorm(log(mean_obs_L)-sqr(sigma_kappa)/2.0, sigma_kappa));
	 // }
	  for (int t=1; t<NY; t++)
	  {
	    Type mean_obs_L = m[t]*S*(S*N[t-1] + (2.0-S)*N[t])/2.0;
	    nll_comp[5] -= dnorm(log(L[t]), log(mean_obs_L)-sqr(sigma_kappa)/2.0, sigma_kappa, true);
	    //SIMULATE{
	     // L[t] = exp(rnorm(log(mean_obs_L)-sqr(sigma_kappa)/2.0, sigma_kappa));
	    //}
	  }
	}
  // Instead of a popcorn model we make it more of a catchability model, I don't think this is a very good model, but we can see what happens
	if(popcorn == 0)
	{
	  for (int t=0; t<NY; t++)
	  {
	  Type mean_obs_m = S*m[t]; 
	  nll_comp[5] -= dnorm(log(mobs[t]), log(mean_obs_m)-sqr(sigma_kappa)/2.0, sigma_kappa, true);
	  //SIMULATE{
	   // mobs[t] = exp(rnorm(log(mean_obs_m)-sqr(sigma_kappa)/2.0, sigma_kappa));
	  //}
	  }
	}
	
	// effort dynamic eq for catch, nothing for C[0]
//Passing in "enable_catcheq" as an indicator: 0 means excluding equation(9), 1 means including equation(9).	
  if(enable_catcheq == 1)
  {
    Type denom = a*B[0]/2.0;
    for (int t=1; t<NY; t++)
    {
      Type mean_obs_C = C[t-1]/B[t-1]*B[t]*pow(B[t-1]/denom,chi);
      nll_comp[6] -= dnorm(log(C[t]), log(mean_obs_C)-sqr(sigma_C)/2.0, sigma_C, true);
      //SIMULATE{
      //  C[t] = exp(rnorm(log(mean_obs_C)-sqr(sigma_C)/2.0, sigma_C));
      //}
    }
  }

	//----------------------------------------------------------------------------
	// Outputs
	//----------------------------------------------------------------------------
	
	REPORT(nll_comp);
  //SIMULATE{
  //REPORT(log_R);
  //REPORT(log_m);
  //REPORT(log_B);
  //REPORT(mobs);
  //REPORT(C);
  //REPORT(L);
  //REPORT(I);
  //REPORT(IR);
  //}
  
	ADREPORT(sigma_tau);
	ADREPORT(sigma_phi);
	ADREPORT(sigma_m);
	ADREPORT(sigma_epsilon);
	ADREPORT(sigma_upsilon);
	ADREPORT(sigma_kappa);
	ADREPORT(sigma_C);
	ADREPORT(q_I);
	ADREPORT(q_R);
	ADREPORT(S);
	ADREPORT(a);
	ADREPORT(chi);
	
	ADREPORT(B);
	ADREPORT(R);
	ADREPORT(m);
	
	Type nll = nll_comp.sum();
	return nll;

}
