clear all; clc; close all;
myPath = ['/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L6'...
    ' - Ambiguity resolution and atmospheric corrections/'];
cd(myPath); addpath(genpath(myPath));
ion.alpha = [2.0724D-08  1.4931D-08 -2.6513D-07 -4.5973D-07];
ion.beta = [1.3581D+05 -1.2248D+05 -1.3794D+06 -2.1847D+06];
Time = 90000;
IP = [ion.alpha; ion.beta];
D_Iono = Klobuchar(Time, IP, lat, long, Az, El);
% Time - seconds of week
% IP(2,4) - alfa and beta-parameters in [sec sec/SC sec/SC^2 sec/SC^3]
% lat - latitude of receiver in radians
% long - longitude of receiver in radians
% Az - azimuth to satellite in radians
% El - elevation angle in radians
% OUTPUT:
% D_Iono - Delay in metres
laloh = fGC2GL(XYZ);
% INPUT:
% XYZ - Cartesian coordinates 
%
% OUTPUT: 
% laloh - [Latitude (rad) Longitude (rad) Height (m)] on ellipsoid WGS 84
[A,E] = fAzimElev(satEcef,recEcef);
% INPUT :
% satEcef - satellite coordinates, earth centred earth framed (vector)
% recEcef - receiver coordinates, earth centred earth framed (vector)
%
% OUTPUT: 
% A - Azimuth angle in radians  
% E - elevation angle in radians