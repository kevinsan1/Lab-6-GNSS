%% covariance calculation
clear all;
X=[ 4.0 	4.2 	3.9 	4.3 	4.1;...
    2.0 	2.1 	2.0 	2.1 	2.2;...
    0.60 	0.59 	0.58 	0.62 	0.63];
S=cov(X')
[afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(X,S,6)