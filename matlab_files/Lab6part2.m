clear all;clc
addpath(genpath(['/Users/kevin/SkyDrive/KTH Work/Period 4',...
    ' 2014/GNSS/Labs/L6 - ',...
    'Ambiguity resolution and atmospheric corrections/']))
load lab4data.mat;
%% Clear unneeded variables
clearvars -except X weightMatrix cofactor
%% Define ambiguities and covariance
an = 4:8;
qn = 4:8;
a = X(an); % float ambiguities from Lab 4
Q = cofactor(qn,qn); % Covariance
clear qn an
%% Use LAMBDA
method = 1;
[afixed,sqnorm,Ps,Qzhat,Z,nfixed,mu]=LAMBDA(a,Q,method);
%% Perform Ratio Test
ratiotest = sqnorm(1)/sqnorm(2);
if ratiotest <= mu
    fprintf('Ratio Test Passed \n')
    fprintf('%0.4f <= %0.4f \n',ratiotest,mu)
else
    fprintf('Ratio Test Failed \n')
    fprintf('%0.4f > %0.4f\n',ratiotest,mu)
end
%% Demo
% [a_ILS,sqnorm]                    = LAMBDA(a,Q,1);
% [a_ILS2]                          = LAMBDA(a,Q,2);
% [a_R]                             = LAMBDA(a,Q,3); 
% [a_B]                             = LAMBDA(a,Q,4); 
% [a_PAR,snt,PsPAR,Qzpar,Zpar,nPAR] = LAMBDA(a,Q,5); 
% [a_RT,snt,Ps,Qzhat,Z,nRT,mu]      = LAMBDA(a,Q,6); 

