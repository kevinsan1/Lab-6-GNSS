clear all;close all;clc;
load troposphereAndIonosphereData.mat
%% Plot troposphere vs angle
%%
averageHeight = sum(height)/length(height);
elevationAngle = 0:90;
P=955;
troposphericCorrection = 0.0022768 * ...
    P/(1-0.00266*cos(2*elevationAngle')...
    -2.8*10^-7*averageHeight);
figure(1);clf;
plot(elevationAngle,troposphericCorrection(:),'.')
%% Plot ionosphere vs angle
plot(el*180/pi,corrIonosphere,'.')
%%
% ionalpha = [2.0724D-08  1.4931D-08 -2.6513D-07 -4.5973D-07];
% ionbeta = [1.3581D+05 -1.2248D+05 -1.3794D+06 -2.1847D+06];
% time = 2 * 24 * 3600;
% iP = [ionalpha; ionbeta];
% for i = 1:length(lat)
%     ionosphericCorrection(i) = Klobuchar(time, iP, lat(i), long(i), az(i), el(i));
% end
% figure(2);
% plot(el*180/pi,ionosphericCorrection,'.')