%% Computation of receiver's position
clear all;
lab1and2Path = ['/Users/kevin/SkyDrive/KTH Work/Period 4 2014/GNSS/Labs/L6'...
    ' - Ambiguity resolution and atmospheric corrections/'...
    'matlab_files/'];
addpath(genpath(lab1and2Path))
cd(lab1and2Path)
c = 299792458; % speed of light (m/s)
mu = 3.986005e14; % universal gravitational parameter (m/s)^3
omega_e_dot = 7.2921151467e-5; % earth rotation rate (rad/s)
f = -4.442807633e-10; % s/m^1/2
time = [2,1,14,0]; % days, hours, minutes, seconds
junk = num2cell(time);
[nday,nhours,nminutes,nseconds] = junk{:};
%% Import observer file
lov033b = importObserverFileAsString('0lov033b.04o', 1, 5629);
% Import P1 numbers and satellite numbers
[rowInObs,nOfRows] = findTimeInObsFunction( lov033b,time ); % match your time with observer time
p1_numbers = importObsP1numbers('0lov033b.04o', rowInObs+1,rowInObs+nOfRows*2); % Import P1 numbers from the matched time above
satelliteNumbers = importObsSatelliteNumbers('0lov033b.04o', rowInObs,rowInObs); % Doesn't change
[xA0,yA0,zA0] = sampleFunction(lov033b); % Record Approximate Position
approxPos = [xA0,yA0,zA0];

%% Import navigation file
navfiles = importNavigationFiles('0lov033b.04n');
%% Match satellite numbers with available satellites
satNumMatch = navfiles(1:8:96,1); % Order of satellite numbers import
%%
sortedSatelliteNumbers = sortrows([satelliteNumbers',p1_numbers],1);
%% Main loop steps 1-14
% Calculates variables needed for correction iterations
count = 1;
% lmat = ones(length(satNumMatch),1);
% amat = ones(length(satNumMatch),4);
% xs = ones(length(satNumMatch),1);
% ys = ones(length(satNumMatch),1);
% zs = ones(length(satNumMatch),1);
% rho = ones(length(satNumMatch),1);
% p1 = ones(length(satNumMatch),1);
% dtsL1_with_dtr = ones(length(satNumMatch),1);
% tAtoS = ones(length(satNumMatch),1);
% v = ones(length(satNumMatch),2)
%%
for i = 1:length(satNumMatch)
    %% Steps 1-14 done inside
    if cell2mat(satNumMatch(i))==sortedSatelliteNumbers(count,1)
        [ lmat(count,1), ...
            amat(count,1:4),...
            rho(count,1),...
            xs(count,1),ys(count,1),zs(count,1),...
            p1(count,1),...
            dtsL1_with_dtr(count,1),...
            tAtoS(count,1),...
            phik(count,1),omegak(count,1),corrIonosphere(count,1),...
            corrTroposphere(count,1),az(count,1),el(count,1),...
            height(count,1),lat(count,1),long(count,1)]... % addition for lab6
            = satLandP( i,sortedSatelliteNumbers(count,2),navfiles,...
            xA0,yA0,zA0,nday,nhours,nminutes,nseconds);
        count = count + 1;
    else
        fprintf('No data for Satellite%3d\n',cell2mat(satNumMatch(i)))
    end
end

%% 17. Repeat steps 11 -16 until the solution has converged.
% The solution has converged if the condition is fulfilled
for i = 1:10
    changeX = (amat'*amat)\(amat'*lmat); % eq. (21)
    v(:,i) = -amat*changeX + lmat; % eq. (17)
    newXYZ = [xA0,yA0,zA0] + changeX(1:3)'; % eq. (22)...
    %     estimated coordinates
    newxyzcell = num2cell(newXYZ);
    [xA0,yA0,zA0] = newxyzcell{:};
    clear newxyzcell;
    rho = sqrt(... % recompute rho
        (xs - xA0 + omega_e_dot * yA0 * tAtoS).^2 + ... % x^2
        (ys - yA0 - omega_e_dot * xA0 * tAtoS).^2 + ... % y^2
        (zs - zA0).^2   ... % z^2
        );
    amat = [-(xs - xA0)./rho,... % recompute matrix A
        -(ys-yA0)./rho,...
        -(zs-zA0)./rho,...
        rho./rho];
    lmat = p1 - rho + c*dtsL1_with_dtr...
        - corrIonosphere - corrTroposphere; % recompute matrix L
    %% 17. Convergence condition
    if i>1 % check for convergence condition
        condition = abs(v(:,end)'*v(:,end)...
            -v(:,end-1)'*v(:,end-1)); % condition must be 1e-5 or lower
        if condition < 1e-4
            fprintf('Convergence condition met = %d\n',condition);
            %% Lab 2
            % Find sigma and Q
            q = inv((amat'*amat));
            sigma_0 = sqrt(v(:,end)'*v(:,end)/(length(amat)-length(q)));
            sigma_x = sigma_0*sqrt(q(1,1));
            sigma_y = sigma_0*sqrt(q(2,2));
            sigma_z = sigma_0*sqrt(q(3,3));
            sigma_t = sigma_0*sqrt(q(4,4))/c;
            pDOP = sqrt(q(1,1)^2 + q(2,2)^2 + q(3,3)^2);
            fprintf('X = %7.3f, mX = %7.3f\n',xA0,sigma_x);
            fprintf('Y = %7.3f, mX = %7.3f\n',yA0,sigma_y);
            fprintf('Z = %7.3f, mX = %7.3f\n',zA0,sigma_z);
            fprintf('T = %0.10f, mt = %0.10d\n',-changeX(4)'/c,sigma_t);
            fprintf('xA0 = %7.3f\n',approxPos(1));
            fprintf('yA0 = %7.3f\n',approxPos(2));
            fprintf('zA0 = %7.3f\n',approxPos(3));
            fprintf('PDOP = %7.3f\n',pDOP);
            fprintf('PDOP = %7.3f\n',...
                sqrt(sigma_x^2 + sigma_y^2 + sigma_z^2));
            break;
        end
    end
end
save troposphereAndIonosphereData.mat;