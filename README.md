# SMC2 sampled mixtures of Gaussian process experts
Implementation of SMC2 mixtures of Gaussian process experts. This is a repository for the software used to create some of the results for the paper 
"Mixtures of Gaussian Process Experts with SMC2'' ([arxiv.org/abs/2208.12830](https://arxiv.org/abs/2208.12830)).

The software approximates posterior distributions for parameters of mixtures of Gaussian process experts. For a detailed description of the sampling procedure see the aforementioned paper.

The results are computed in three steps. First, SMC2 is used to sample posterior distributions of the model parameters using the main scripts found in the root folder. Second, densities scripts, again in the project root folder, are used to
compute the predictive distributions, posterior similarity matrices, and posterior for the number of used experts. These results can be visualized with the respective plotting scripts. For the higher-dimensional exaples, the posterior similarity matrices
need to be reorganized with 
the supplied R script.

# Installation
Clone or download the repository and add the include folder and subfolders to your MATLAB path. This should happen automatically on running the main script in the project root folder.

# References
If you find the software useful, please cite [arxiv.org/abs/2208.12830](https://arxiv.org/abs/2208.12830).
