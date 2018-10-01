# abc-toy-problems
A series of 4 toy problems to demonstrate the utility and function of likelihood-free Bayesian inference, Approximate Bayesian Computation (ABC). These problems were created as part of my Master of Research (MRes) thesis on 'Approximate Bayesian Tomography', applying the ideas of likelihood free inference to geophysics.

Toy Problem 1: Estimates the parameters to a Gaussian distribution (mean and standard devation) using an ABC rejection sampler

Toy Problem 2: Performs linear regression using ABC-MCMC

Toy Problem 3: Estimates the parameters to a bivariate Gaussian distribution using ABC-MCMC

Toy Problem 4: Estimates the parameter to a 'banana distribution' (a twister bivariate Gaussian distribution) by introducing Adaptive Metropolis and Delayed Rejection to ABC-MCMC

The sampling algorithms are found in /dram. This code is all built on top of an [MCMC toolbox by Marko laine](http://helios.fmi.fi/~lainema/mcmc/).
