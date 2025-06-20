# advection-diffusion-neg

This repository provides code for 
- simulating the evolution equation and
- visualizing the results of eigenvalue analysis

in the following paper. I would appreciate your citing the following paper when you publish your results using this code.

Ohtake, K. (2024). Pattern formation by advection-diffusion in new economic geography. arXiv preprint arXiv:2407.05804.  
<a href="https://doi.org/10.48550/arXiv.2407.05804" target="_blank" rel="noopener noreferrer">https://doi.org/10.48550/arXiv.2407.05804</a>

# Overview

## advection_diffusion_model_autoparam.ipynb
Language: 
Julia ver 1.10.2  
Packages:  
CSV ver. 0.10.13  
DataFrames ver. 1.6.1  
Distributions ver. 0.25.107  
Format ver. 1.3.7  
IJulia ver. 1.24.2  
Plots ver. 1.40.2  

### How to use

In the second cell, specify the values of $\tau>0$ as the list used in the outer loop.  
Then, specify parameters $\mu\in(0,1)$ and $\sigma>1$ in the list used in the inner loop.  
For a set of parameters, five simulations are performed for randomly generated initial values.  

## eigenvalue_QLLU-AD/heatmap_sigma.py
This is a code for computing eigenvalues of QLLU-AD and drawing heat maps.  
Similarly, other directories named "eigenvalue_X" contain Python programs that plots heatmaps for X model.

Language: Python ver 3.12.4  
Packages:  
matplotlib ver 3.10.1  
numpy 2.2.3  

## eigenvalue_QLLU-AD/contour_sigma.py

This is a code for computing eigenvalues of QLLU-AD and plotting critical curves.  
Similarly, other directories named "eigenvalue_X" contain Python programs that plots critical curves for X model.

Language: Python ver 3.12.4  
Packages:  
matplotlib ver 3.10.1  
numpy 2.2.3
