# MFG

The approach of mean-field games (MFG) were introduced by J.-M. Lasry and P.-L. Lions and were successfully adapted to the economic problems.
Mathematically, mean field equilibrium (i.e., Nash point for an infinite number of players) leads to the system of two parabolic partial differential equations: Kolmogorov and Hamilton-Jacobi-Bellman ones coupled by some differential equality.

This implemnetaion of numerical solving of Kolmogorov equaiton for 1D case.

![Eq1](https://latex.codecogs.com/gif.latex?%5Cfrac%7B%5Cpartial%20m%7D%7B%5Cpartial%20t%7D%20-%20%5Cfrac%7B%5Csigma%5E2%7D%7B2%7D%5Cfrac%7B%5Cpartial%5E2%20m%7D%7B%5Cpartial%20x%5E2%7D%20&plus;%20%5Cfrac%7B%5Cpartial%20%28%5Calpha%20m%29%7D%7B%5Cpartial%20x%7D%20%3D%200%5C%5C%5C%5C%20m%28x%29%3Dm_0%28x%29%5Cquad%5Cforall%5C%2Cx%20%5Cin%20%280%2C1%29%5C%5C%5C%5C%20%5Cfrac%7B%5Cpartial%20m%7D%7B%5Cpartial%20x%7D%28t%2C0%29%3D%5Cfrac%7B%5Cpartial%20m%7D%7B%5Cpartial%20x%7D%28t%2C1%29%3D0%20%5Cquad%20%5Cforall%5C%2Ct%20%5Cin%20%280%2CT%29%5C%5C%5C%5C%20%5Calpha%28t%2C0%29%3D%5Calpha%28t%2C1%29%3D0%20%5Cquad%20%5Cforall%5C%2Ct%20%5Cin%20%5B0%2C%20T%5D)
