Workflow: 

1. Sampler.C allows the user to pass in parameters and initialize bivariate Gaussian copula model with Elliptic Power or Gaussian marginal distributions (for v_n and v'_n). It also uses ROOT's built in monte carlo sampling method to obtain vectors for v_n and v'_n by sampling the distribution. 
 
2. Quantities.C allows the user to evaluate >30 observables using the vectors with v_n and v'_n.

3. Datapoints.C allows the user to fix certain parameters, and then iterates through a set of parameter values for rho (the pearson's correlation coefficient between v_n and v'_n) and the parameter governing the width of the v'_n distribution. At each point in this 2d phase space, it generates a distribution with the 5 parameters, and calls Sampler.C to evaluate random datapoints from that distribution. Then Quantities.C is called, and evaluates each observable using the sampled data. The values for each of these quantities is stored in a Root file.

4. Dataread opens the root file, and uses the subsampling method to determine errors, by reading each datapoint from the root file, and determining for each observable, the mean and standard error. This is written to a more compact root file that can be used for plotting data etc. 
