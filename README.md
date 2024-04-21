Repository for the paper "Computation of two-dimensional Stokes flows via lightning and AAA rational approximation" by <a href="https://yidanxue.github.io">Yidan Xue</a>, <a href="https://www.maths.ox.ac.uk/people/sarah.waters">Sarah L. Waters</a> and <a href="https://people.maths.ox.ac.uk/trefethen/">Lloyd N. Trefethen</a>. The arXiv preprint can be found at https://arxiv.org/abs/2306.13545. The SISC publication can be found at https://doi.org/10.1137/23M1576876.

The codes need to be run in MATLAB or Octave, its free alternative. The recommended version is MATLAB_R2023b. The user needs to have aaa.m in the same folder (or add it in the default MATLAB path) for the AAA rational approximation. The easist way to do this is to download Chebfun at https://www.chebfun.org. 

There are 4 MATLAB codes in this repository to reproduce the results presented in our paper:
1) 'constricted_channel.m' computes 2D Stokes flows through a channel with a smooth constriction using the AAA algorithm. One can change the lambda value to see the changes in poles locations and flow characteristics for different constriction amplitude (or expansion amplitude for negative lambda). One can also turn off 'AAA' by changing lines 35-39 to 'Pol={};' to see how that will change the convergence of the computation. This code reproduces Fig. 3 in our paper.
2) 'Stokes_flow_between_two_cylinders.m' computes 2D Stokes flows between two rotating and translating cylinders using a series method. This code presents how to compute Stokes flows in a multiply-connected domain using a series method. By changing parameter values based on Table 1 (the current setting is for Fig. 6a), one can reproduce Fig. 6 in our paper. This code can also reproduce Fig.8 (ellipse in ellipse) with minor changes.
3) 'particle_in_bifurcation.m' reproduces Fig. 11 in our paper. This code shows how to combine the lightning algorithm, the AAA algorithm and a series method to compute 2D Stokes flows in general domains. One can adapt this code for their own Stokes computations.
4) 'heart_eddies.m' reproduces Fig. 10. This code demonstrates the accuracy of our computation for visualising Moffatt eddies at 13-digit.

In the tutorial folder, there are additional tutorials to work you through lightning, AAA, and series methods for rational approximation of a variety of Stokes problems. The best way to compute these tutorials is to use the <tt>publish</tt> function in MATLAB.

We will keep adding new MATLAB codes in this repository for more 2D Stokes flow computations.
