function Y = EOM_func(t,X)
%EOM_FUNC 
% This function serves at the differential equations to be integrated
% elsewhere in the folder to compute the interplanetary trajectory
% 
% Inputs **
% t = Scenario time in julian days
% X = state vector
%       [x, y, z, xDot, yDot, zDot, centralPlanet]
%
% Outputs **
% Y = state vector to be integrated

% All values are in km^3/sec^2
GM1 = 22031.868551; % Mercury GM
GM2 = 324858.592000; % Venus GM
GM3 = 398600.435507; % Earth GM
GM4 = 42828.375816; % Mars GM
GM5 = 126712764.100000; % Jupiter GM
GM6 = 37940584.841800; % Saturn GM
GM7 = 5794556.400000; % Uranus GM
GM8 = 6836527.100580; % Neptune GM
GM9 = 975.500000; % Pluto GM
GMS = 132712440041.279419; % Sun GM

mu = [GM1, GM2, GM3, GM4, GM5, GM6, GM7, GM8, GM9, GMS];


R = X(1:3);
RmagSat = norm(R);
Y = 0.*X;
% Velocity
Y(1) = X(4);
Y(2) = X(5);
Y(3) = X(6);

% Acceleration
Y(4) = -mu(X(7))/RmagSat^3 * X(1);
Y(5) = -mu(X(7))/RmagSat^3 * X(2);
Y(6) = -mu(X(7))/RmagSat^3 * X(3);

% if X(7) == 10
%     % (N-1) Body Perturbations
%     for i = 1:10
%         if i == X(7)
%             continue
%         else 
%             t0 = datetime(t/3600/24,'ConvertFrom','juliandate', ...
%                 'Format',"MMM dd, yyyy HH:mm:ss.SSS");
%             ephTime = cspice_str2et(char(t0));
%             [planetRV, ~] = cspice_spkezr(char(string(i)),ephTime,'J2000','NONE','10');
%             Rplanet = (planetRV(1:3));
%             Rmag = norm(Rplanet);
%             Runitvec = Rplanet/Rmag;
% 
%             Y(4) = Y(4) + mu(i)/Rmag^3*(3*dot(Runitvec,R)*Runitvec(1) - R(1));
%             Y(5) = Y(5) + mu(i)/Rmag^3*(3*dot(Runitvec,R)*Runitvec(2) - R(2));
%             Y(6) = Y(6) + mu(i)/Rmag^3*(3*dot(Runitvec,R)*Runitvec(2) - R(3));
%         end
%     end
% end
