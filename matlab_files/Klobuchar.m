function D_Iono = Klobuchar(Time, IP, lat, long, Az, El )
% INPUT: 
% Time - seconds of week
% IP(2,4) - alpha and beta-parameters in [sec sec/SC sec/SC^2 sec/SC^3]
% lat - latitude of receiver in radians
% long - longitude of receiver in radians
% Az - azimuth to satellite in radians
% El - elevation angle in radians
% OUTPUT:
% D_Iono - Delay in metres
%
% All input angles should be in radians
% help Klobuchar;
% REFERENCES: Hoffman-Wellenhof (2001) approach page 103
%             Klobuchar (1987)
%             A Leick (2004)
%
% Last updated by MH 2014-04-24

% Constants *****************************************************
    C=299792458; % m/s

% INPUT:-----------------------------------------------------------------
    if IP == 0
        D_Iono = 0;
        return;
    end

    alpha = IP(1,:);
    beta = IP(2,:);

    El = El/pi; % semicircles; % JV 2005-04-14
    %Az = Az; % radians

%****************************************************************
% Table 6.3 The Broadcast Ionospheric Model and Klobuchar 1987

% Calculate the earth centred angle:
    Psi = 0.0137/(El+.11) - 0.022; % semicircles
    
% Compute the subionospheric latitude: Ionospheric Point (IP)
    
    phiIP = lat+Psi*cos(Az); % semicircles
   
    if phiIP > 0.416;          % semicircles
        phiIP = 0.416;         % semicircles
    elseif phiIP < -0.416;     % semicircles
        phiIP = -0.416;        % semicircles
    end

% Compute the subionospheric longitude:
    lambdaIP = long+(Psi*sin(Az))/cos(phiIP*pi);% semicircles
          
% Compute the geomagnetic latitude
    cGeoMagneticLatIP = phiIP+0.064*cos((lambdaIP-1.617)*pi);% semicircles
    
% Find the local time
    t = lambdaIP*43200 + Time;   
    t = mod(t,86400);
    

% Amplitude of correction curve
    AMP = alpha(1) + alpha(2)*cGeoMagneticLatIP + alpha(3)...
        *cGeoMagneticLatIP^2 + alpha(4)*cGeoMagneticLatIP^3;
    if AMP<0
        AMP = 0;
    end

% Period of the correction curve
    PER =  beta(1) + beta(2)*cGeoMagneticLatIP...
        + beta(3)*cGeoMagneticLatIP^2 + beta(4)*cGeoMagneticLatIP^3;
    if PER<72000,   PER = 72000;  end

% Phase of correction curve
    x = 2*pi*(t-50400)/PER;
   
% Compute the actual slant factor
    F = 1+16*(0.53-El)^3; 

% Then compute the Ionospheric time delay
      
    % According to ICD-GPS-200C pp. 125 - 127
        if abs(x)>=1.57
            DT_Iono=F*5e-9;
        else
            DT_Iono=F*(5e-9 + AMP*(1-x^2/2+x^4/24));
        end

% Transform the time delay to distance delay
    D_Iono = DT_Iono*C;