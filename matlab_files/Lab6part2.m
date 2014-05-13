clear all;clc
load lab4data.mat;
%% 
ahat = X(4:end); % float ambiguities from Lab 4
dahat = [ahat,L(1:5)]';
m=5;
method = 3; % Integer rounding method
xhat = mean(dahat,2);
%%
Xc = dahat-xhat*ones(1,m); % centered data matrix
%%
Qahat = (1/m)*Xc'*Xc; % covariance matrix
%%
%%
[afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(ahat,Qahat,method)