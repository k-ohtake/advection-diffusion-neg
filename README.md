# advection-diffusion-neg

This repository provides simulation codes for an evolution equation discussed in the following paper.  
I would appreciate your citing the following paper when you publish your results using this code.

Ohtake, K. (2024). Pattern formation by advection-diffusion in new economic geography. arXiv preprint arXiv:2407.05804.  
<a href="https://doi.org/10.48550/arXiv.2407.05804" target="_blank" rel="noopener noreferrer"> https://doi.org/10.48550/arXiv.2407.05804</a>

# Overview

Language: 
Julia ver 1.10.2  
Packages:  
CSV ver. 0.10.13  
DataFrames ver. 1.6.1  
Distributions ver. 0.25.107  
Format ver. 1.3.7  
IJulia ver. 1.24.2  
Plots ver. 1.40.2  

# How to use

In the second cell, specify the values of $\tau>0$ as the list used in the outer loop.  
Then, specify parameters $\mu\in(0,1)$ and $\sigma>1$ in the list used in the inner loop.  
For a set of parameters, five simulations are performed for randomly generated initial values.  
