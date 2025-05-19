function C = mgaPlotterV2(mission,planetStart,rDepart,C3,planetTarget, ...
                          targetOrbitShape,maxN,X)
% mgaPlotterV3 - a derivative of mgaDesignV2, the fitness function used for
% a genetic algorithm, designed to plot the traejctories as identified by
% the genetic algorithm
%
% Inputs**
% mission = Uranus/Galileo. Effects the types of plot to output
% planetStart = planet at the start of the trajectory
% rDepart = paarking orbit radius
% C3 = injection energy from inertial upper stage
% planetTarget = target destination
% targetOrbitShape = desired orbit upon arrival to target planet
% maxN = maximum gravity assists
% X(1) = number of gravity assists
% X(2:maxN+1) = gravity assist planets
% X(maxN+2) = departure date (julian days)
% X(maxN+3:2*maxN+2) = times of flight between gravity assists (days)
% X(2*maxN+3) = final leg time of flight (days)
% X(2*maxN+4) = angle 1 for possible resonant orbit (radians)
% X(2*maxN+5) = angle 2 for possible resonant orbit (radians)
%
% Outputs**
% C = cost function. Sums total deltaV cost and adds penalties for low
% energy transfers or periapses that fly within their planet
%
%
% Planet numbers and the number of flybys should be integers
%
% References: [13-15, 17, 18, 21, 22, 27]


% Variables to be optimized
t0 = X(1); % launch time
numGA = X(2); % number of gravity assits
gaPlanets = X(3:2+numGA); % gravity assist order
TOFs = X((maxN+3):(2+maxN+numGA)); % times of flight for each leg
tf = X(2*maxN+3); % time for the final leg of the trajectory
ang = X((2*maxN+4):(2*maxN+5)); % angle of fly-by during resonance

% Load JPL HORIZON kernals
cspice_furnsh('de440.bsp');
cspice_furnsh('s970311a.bsp');
cspice_furnsh('s980326a.bsp');
cspice_furnsh('latest_leapseconds.tls');

% Planetary Mass Parameters [17, 18]
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

% Planetary Distances [17, 18]
au = 1.495978e8; % km, 1 au
mercuryAU = 0.3871; % au
venusAU = 0.7233; % au
earthAU = 1; % au
marsAU = 1.5237; % au
jupiterAU = 5.2028; % au
saturnAU = 9.5388; % au
uranusAU = 19.1914; % au
neptuneAU = 30.0611; % au
plutoAU = 39.5294; % au
planetDistances = [mercuryAU, venusAU, earthAU, marsAU, jupiterAU, ...
                   saturnAU, uranusAU, neptuneAU, plutoAU].*au;

% Compute the radii of the spheres of influence
for l = 1:length(mu)-1
    rSoi_planets(l) = (mu(l)/mu(10))^(2/5)*planetDistances(l);
end

% Planetary Mean Radii
rMercury = 2439; % km
rVenus = 6051.8; % km
rEarth = 6371; % km
rMars = 3389.5; % km
rJupiter = 69911; % km
rSaturn = 58232; % km
rUranus = 25362; % km
rNeptune = 24622; % km
rPluto = 1188.3; % km
rPlanet = [rMercury, rVenus, rEarth, rMars, rJupiter, ...
           rSaturn, rUranus, rNeptune, rPluto];

% List of the planets to be visited in this trajectory
planetVisited = [planetStart gaPlanets(1:numGA) planetTarget];

% List of the durations between planets
legDurations = [0 TOFs tf];

% List of the times in julian days
planetEpochs = t0.*ones(length(planetVisited),1);
for i = 1:numGA
    planetEpochs(i+1:end) = planetEpochs(i+1:end) + TOFs(i);
end

% For the terminal leg, add the final duration
planetEpochs(end) = planetEpochs(end) + tf;

% Create an N by 1 cell, enough to containt the states of the planets
% visited at the times assigned by the genetic algorithm
planetStates = cell(length(planetVisited),1);

% For each planet visited, perform the following...
for i = 1:length(planetVisited)
    % Convert the date from julian to a format recognized by the SPICE
    % toolkit
    tDate = datetime(planetEpochs(i),'ConvertFrom','juliandate',...
        'Format',"MMM dd, yyyy HH:mm:ss.SSS");
    ephTime = cspice_str2et(char(tDate));
    
    % Get the state of the planets visited
    [state, ~] = cspice_spkezr(char(string(planetVisited(i))), ...
        ephTime,'J2000','NONE','10');

    planetStates{i} = state;
end

% In case of resonance orbits, this is the counter. There should never be
% more than 2 resonance orbits per trajectory
resonanceCounter = 0;

% Create a N by 2 cell to store the inbound and outbound velocities during
% each gravity assist
lambertSolutions = cell(length(planetVisited),2);

% For every planet visited, except the target planet, perform the following
for i = 1:length(planetVisited)-1
    % If the current planet is the same as the next planet, check to see if
    % its trajevtory is a resonance orbit
    if planetVisited(i) == planetVisited(i+1)
        % Get the mean time for the current planet to orbit the sun in days
        Tplanet = sqrt(4*pi^2*planetDistances(planetVisited(i))^3 ...
            /mu(10))/3600/24; % days
        
        % The possible resonance orbits being considered are 1:2, 1:3, 2:2,
        % and 2:3. If the orbit does not fall within +/- 2 days of these
        % ratios, then the orbit is not in resonance

        % Number of satellite orbits
        satRotations = [1 2];
        % Number if planet orbits
        planetRotations = [2 3];

        % check each combination of satellite roations and planet rotations
        for j = satRotations
            for k = planetRotations
                % Difference in leg TOF and a resonance orbit
                resonantDeltaT = abs(legDurations(i+1) - k/j*Tplanet);
    
                % If the absolute difference in days is less than 2,
                % and this is not the third or more resonance orbit
                % this trajectory
                if resonantDeltaT <= 2 && resonanceCounter <= 2
                    % Get the number of satellite orbits
                    if j == k
                        % If the satellite revolutions equals the planet
                        % revolutions, then only one orbit has occured
                        revs = 0;
                    else
                        revs = j;
                    end

                    % Add to the resoannce orbit counter
                    resonanceCounter = resonanceCounter + 1;
                    break
                else
                     % Otherwise, the trajectory is not close enough to be 
                    % resonance and casue cause trouble for a lambert's 
                    % problem solver
                    revs = 0;
                    if resonanceCounter == 3
                        % If the issue is that this is the third resonant
                        % orbit, set the score to a large value and return
                        % to the function call within the genetical
                        C = 1000;
                        return
                    end
                end
            end
            % If the inner-loop was broken, break the outer-loop
            % becuase the number of orbits was determined
            if revs ~= 0
                break
            end
        end               
    else
        revs = 0;
    end

    % If the next leg of the trajectory is not resonant, then solve
    % Lambert's problem as normal and save the departing velocity and the
    % arrival velocity, from their repsective planets
    if revs == 0
        [vOut, vIn] = glambert(mu(10), ...
                               planetStates{i}, ...
                               planetStates{i+1}, ...
                               legDurations(i+1)*24*3600, 0);

        lambertSolutions{i,2} = vOut;
        lambertSolutions{i+1,1} = vIn;

    % Otherwise, perform the following...
    else
        % If currently looking at the departing planet to the first gravity
        % assist, try to compute the solutions to lambert's problem with
        % the allocated number of complete revolutions
        if i == 1
            try
                [vOut, vIn] = glambert(mu(10), ...
                                   planetStates{i}, ...
                                   planetStates{i+1}, ...
                                   legDurations(i+1)*24*3600, revs);
            
                % If there is no solution possible, recompute but set the
            % revolutions to 0. Because this is a powered maneuver, the
            % first orbit should not be resonant but this is just to cover
            % all possible scenarios
            catch
                [vOut, vIn] = glambert(mu(10), ...
                                   planetStates{i}, ...
                                   planetStates{i+1}, ...
                                   legDurations(i+1)*24*3600, 0);
            end

            % Save the inital velocity in the current row, second column
            % while the terminal velocity in the first column, next row
            lambertSolutions{i,2} = vOut;
            lambertSolutions{i+1,1} = vIn;
        else
            % If not on the first planet, try to compute the solutions to
            % lambert's problem
            try
                % This method should not fail, but if it does save empty
                % arrays in the corresponding location
                [vOut, vIn] = glambert(mu(10), ...
                               planetStates{i}, ...
                               planetStates{i+1}, ...
                               legDurations(i+1)*24*3600, revs);
                lambertSolutions{i,2} = vOut;
                lambertSolutions{i+1,1} = vIn;
            catch
                % There is additional code to take care of trajectory legs
                % that are resonant and the velocities between resonant
                % orbit passes are unknown
                lambertSolutions{i,2} = [];
                lambertSolutions{i+1,1} = [];
            end
        end
    end
end

% Tolerance limit to solve for the orbit eccentricities
tol = 1e-9;
% Empty array to store all of the Delta-Vs from Earth departure to orbit
% insertion
deltaVs = zeros(length(planetVisited),1);
% A corresponding array of penalties for poor trajectories
penalties = zeros(length(planetVisited),1);

% For each planet visited except the target planet, perform the following
for i = 1:length(lambertSolutions)-1
    % Get the current planet
    planetNum = planetVisited(i);
    % Get the planet's position and velocitiy states
    planetR = planetStates{i}(1:3);
    vPlanet = planetStates{i}(4:6);
    
    % If the second position of this row is empty, then I know there is a
    % resonant orbit in the next leg
    if isempty(lambertSolutions{i,2})
        % This segment of code should not occur frequently, only when there
        % is a failure by the Lambert's solver to produce a solution. This
        % process is from [22] and it details how to calculate the outbound
        % and inbound velocities between the resonant orbit planet passes

        % Get the heliocentric semi-major axis
        a = ((legDurations(i+1)*24*3600)^2/4/pi^2*mu(10))^(1/3);

        % Get the spacecraft's 1st outbound velocity magnitude
        Vout1 = sqrt(mu(10)*(2/norm(planetR) - 1/a));

        % Depending on the current planet, collect the appropriate relative
        % outbound velcoity
        if i == 1
            vInfIn1 = lambertSolutions{i,2} - vPlanet;
        else
            vInfIn1 = lambertSolutions{i,1} - vPlanet;
        end
        
        % One of two b-plane angles, the second is actually is given by the
        % genetic algorithm
        theta = acos((Vout1^2 - norm(vInfIn1)^2 - norm(vPlanet)^2)/ ...
            (-2*norm(vInfIn1)*norm(vPlanet)));

        % Relative, outbound velocity vector after the first gravity assist
        vInfOut1 = norm(vInfIn1).*[cos(pi-theta);
                            sin(pi-theta)*cos(ang(resonanceCounter));
                           -sin(pi-theta)*sin(ang(resonanceCounter))];

        % Saving the outbound velocity
        lambertSolutions{i,2} = vInfOut1 + vPlanet;

        % Rotation matrix to convert from the velocity-normal-conormal 
        % (VNC) reference frame to the ecliptic plane
        Vhat = vPlanet./norm(vPlanet);
        Nhat = cross(planetR,vPlanet)./norm(cross(planetR,vPlanet));
        Chat = cross(Vhat,Nhat);

        That = [Vhat Nhat Chat]';

        % Relative, outbound velocity in the eccliptic plane
        VInfOutEcc = That*vInfOut1;

        % Relative, inbound velocity in the eccliptic plane
        VInfInEcc = VInfOutEcc + vPlanet - planetStates{i+1}(4:6);

        % Relative, inbound velocity in the VNC frame
        vInfIn2 = That\VInfInEcc;

        % Save the inbound velocity before the second gravity assist of the
        % resonance orbit
        lambertSolutions{i+1,1} = vInfIn2 + planetStates{i+1}(4:6);

        % Old code to compute the trajectory penalities for a low flyby
        % radius
        % aFlyBy = -mu(planetVisited(i))/norm(vInfIn1);
        % eFlyby = eOutCalc(planetNum,vInfIn1,vInfOut1,tol);
        % rP = aFlyBy*(1-eFlyby);
        % 
        % flybyCoeff = 1.05;
        % if rP < flybyCoeff*rPlanet(planetNum)
        %     penalties(i) = -10*log10(rP/(flybyCoeff*rPlanet(planetNum)));
        % else
        %     penalties(i) = 0;
        % end
        % deltaVs(i) = 0;

        % New code to compute the Delta-V across the gravity assist. It
        % should be near zero but if it isn't, then the genetica algorithm
        % will likely ignore this trajectory until the former gravity
        % assist can provide a trajectory with enough energy
        vIn = lambertSolutions{i,1};
        vOut = lambertSolutions{i,2};
        
        [deltaVs(i), penalties(i)] = gadV(...
                                       vIn, vOut, vPlanet, planetNum, tol); 
    
    % If the next orbit is not a resonant orbit, perform the following...
    else
        % If this is the starting planet, perform the following...
        if i == 1
            % Calculate the needed departure excess velocity
            departureVInf = norm(lambertSolutions{1,2} - ...
                planetStates{1}(4:6));
            
            % Compute the excess velocity as a result of the intertial
            % upper stage
            vC3 = sqrt(C3 + 2*mu(planetStart)/rDepart);
            
            % Compute the difference in velocities of what is needed and
            % what can be provided by the inertial upper stage
            deltaVdeparture = vC3 - ...
                        sqrt(departureVInf^2 + 2*mu(planetStart)/rDepart);

            % If the difference is negative, then it means the inertial
            % upper stage did not have enough injection energy and the
            % spacecraft will need to fill the gap
            if deltaVdeparture < 0
                deltaVs(1) = abs(deltaVdeparture);
            
            % If the value is positive, then the inertial upper stage was
            % not used compeltely and no Delta-V is required
            else
                deltaVs(1) = 0;
            end
            

        % If not looking at the departure planet, the perform the
        % following...
        else
            % Get the inbound and outbound velocities for thie gravity
            % assist
            vIn = lambertSolutions{i,1};
            vOut = lambertSolutions{i,2};
            
            % Determine if there needs to be an penalties or maneuvers
            % performed
            [deltaVs(i), penalties(i)] = gadV(...
                vIn, vOut, vPlanet, planetNum, tol);
        end

    end
end

% After computing the data for departure and the intermediate flybys, the
% last step is to calculate the Delta-V needed to enter the provided target
% orbit

% Based on the provided orbit information, get the target orbits periapsis
% and eccentricity
rPTarget = targetOrbitShape(1);
eTarget = targetOrbitShape(2);

% Compute the target orbit's semi-major-axis
aTarget = rPTarget/(1-eTarget);

% Compute the velocity at the periapsis of the target orbit
vPtarget = sqrt((1+eTarget)/(1-eTarget)*mu(planetTarget)/aTarget);

% Compute the inbound relative velocity to the target planet
arrivalVInf = norm(lambertSolutions{end,1} - planetStates{end}(4:6));

% Compute the difference in velocities at the target periapsis
deltaVarrival = abs(vPtarget - ...
                    sqrt(arrivalVInf^2 + 2*mu(planetTarget)/rPTarget));

% Save this difference
deltaVs(end) = deltaVarrival;

% Sum up all Delta-Vs and all penalties to get the score for the trajectory
C = sum(deltaVs) + sum(penalties);

%% Plotting the trajectories identified by the genetic algorithm

% Total number of days traveled
finalTOF = sum(TOFs(1:numGA)) + tf;

% Array of days from t0 to tf
t = t0:1:(t0 + finalTOF);

% Empty array to store galileo's flight data
galileoRV = [];

% Create a list of the unique planets visited on the trajectory
gaPlanets = unique(gaPlanets(gaPlanets > 0));

% empty cell to store all of the positions of the planets from t0 to tf
planetPositions = cell(9,1);

% For each time step in t, get the state of the planets
for i = 1:length(t)
    tDate = datetime(t(i),'ConvertFrom','juliandate','Format', ...
        "MMM dd, yyyy HH:mm:ss.SSS");
    ephTime = cspice_str2et(char(tDate));

    % For each planet, collect the state at the current time
    for j = gaPlanets
        [planetRV, ~] = cspice_spkezr(char(string(j)),ephTime,...
            'J2000','NONE','10');

        planetPositions{j} = [planetPositions{j}; planetRV'];
    end

    % If all gravity assist planets are different then the starting planet,
    % then collect the starting planet's state for time t
    if all(gaPlanets ~= planetStart)
        [planetRV, ~] = cspice_spkezr(char(string(planetStart)),...
            ephTime,'J2000','NONE','10');

        planetPositions{planetStart} = ...
            [planetPositions{planetStart}; planetRV'];
    end
    
    % Collect the state for the target planet
    [planetRV, ~] = cspice_spkezr(char(string(planetTarget)),ephTime,...
        'J2000','NONE','10');
    planetPositions{planetTarget} = ...
        [planetPositions{planetTarget}; planetRV'];

    
end

% Open figure to plot the orbit tracks
figure;
% Color ocdes that closely resemble the colors of each planet
mercuryColor = [122 122 122]./255;
venusColor = [255 109 0]./255;
earthColor = [0 125 235]./255;
marsColor = [161 66 3]./255;
jupiterColor = [227 155 11]./255;
saturnColor = [];
uranusColor = [148 246 247]./255;
neptuneColor = [];
plutoColor = [];

% Cell containing the planet's corresponding colors
planetColors = {mercuryColor venusColor earthColor marsColor ...
    jupiterColor saturnColor uranusColor neptuneColor plutoColor};

% String to containg the planets in the plot
legendSTR = {};

% For each planet, perform the following
for i = 1:9
    % If the cell is not empty, plot the orbit of the planet
    if ~isempty(planetPositions{i})
        plot(planetPositions{i}(:,1),planetPositions{i}(:,2), ...
            'Color',planetColors{i})
        hold on;

        % Add to the dynamic legend depending on the planet
        if i == 1
            legendSTR = [legendSTR "Mercury"];
        elseif i == 2
            legendSTR = [legendSTR "Venus"];
        elseif i == 3
            legendSTR = [legendSTR "Earth"];
        elseif i == 4
            legendSTR = [legendSTR "Mars"];
        elseif i == 5
            legendSTR = [legendSTR "Jupiter"];
        elseif i == 6
            legendSTR = [legendSTR "Saturn"];
        elseif i == 7
            legendSTR = [legendSTR "Uranus"];
        else
            legendSTR = [legendSTR ""];
        end
    end
end

% Dynamic storage for the orbit track of the spacecraft
satT = [];
satRV = [];

% variable to offset the period of time for use in Lambert's problem
deltaT = 0;
% Function options for ode45
option2 = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

% For each gravity assist perform the following...
for i = 1:size(lambertSolutions,1)-1
    % If the current index is not the last...
    if i ~=size(lambertSolutions,1)-1
        % Get a new time t from t0+dektaT to the next leg duration
        t = ((0:1:TOFs(i)) + t0 + deltaT).*24*3600;
    else
        % Otherwise use the terminal leg's time of flight
        t = ((0:1:tf) + t0 + deltaT).*24*3600;
    end

    % Use ode45 to propagate the two-body problem for each leg of the
    % interplanetary trajectory
    [tEnd, Y] = ode45(@EOM_func,t, ...
        [planetStates{i}(1:3);lambertSolutions{i,2}; 10],option2);
    
    % Update the offset
    if i ~= size(lambertSolutions,1)-1
        deltaT = deltaT + TOFs(i);
    end

    % Add to the spacecraft state matrix
    satT = [satT; t'];
    satRV = [satRV; Y];
end

% Plot the spacecraft's trajectory
plot(satRV(:,1),satRV(:,2),'k')

% If the mission isn't to recreate Galileo, perform the folowing
if mission ~= "Galileo"
    legend([legendSTR,"GA Result"],"Location","best")
end
    
axis equal
axis off

% If the mission is to recreate Galileo, perform the following...
if lower(mission) == lower("Galileo")
    % Start of the Galileo mission
    galileoStart = juliandate(datetime(...
        "1989 OCT 19 01:29:33.260",...
        "InputFormat","yyyy MMM dd HH:mm:ss.SSS",...
        "Format","yyyy MMM dd HH:mm:ss.SSS"));
    % Galileo arrival to Jupiter
    galileoEnd =   juliandate(datetime(...
        "1995 DEC 07 21:54:44.782",...
        "InputFormat","yyyy MMM dd HH:mm:ss.SSS",...
        "Format","yyyy MMM dd HH:mm:ss.SSS"));
    
    % If the genetic algorithm's trajectory ends before Galileo on the
    % calendar, then project Galileo's orbit until the end of the
    % spacecraft trajectory
    if galileoEnd < satT(end)/24/3600
        galileoT = galileoStart:1:satT(end)/24/3600;
    else
        % Otherwise, model Galileo from start to finish as normal
        galileoT = galileoStart:1:galileoEnd;
    end

    % For each time step, collect the state of Galileo
    for i = galileoT
        tDate = datetime(i,'ConvertFrom','juliandate','Format', ...
            "MMM dd, yyyy HH:mm:ss.SSS");
        ephTime = cspice_str2et(char(tDate));
        [rv, ~] = cspice_spkezr(char("-77"),ephTime,'J2000','NONE','10');
        galileoRV = [galileoRV rv];
    end
    
    % Plot Galileo's trajectory
    plot(galileoRV(1,:),galileoRV(2,:),'r')
    
    
    legend("Venus","Earth","Jupiter","GA Recreation",...
        "Galileo Flight Data","Location","best")
    
    % Plot of geneitc algorithm's trajectory versus Galileo, position
    % magnitude
    figure
    satMag = vecnorm(satRV(:,1:3)');
    galileoMag = vecnorm(galileoRV(1:3,:));
    plot((satT-satT(1))./24./3600,satMag./au,'k',...
        galileoT-galileoT(1),galileoMag./au,'r');
    xlabel("Time Since Departure (days)")
    ylabel("Distance (AU)")
    title("Magnitude of the Distance From the Sun")
    legend("GA Recreation","Galileo Flight Data","Location","best")
    
    % Plot of genetic algorithm's trajectory versus galileo, differences in
    % position magnitude
    figure
    % If the number of samples in satT is greater than galileoRV, that
    % means there is more data for the spacecraftand so it's trajectory
    % will only be considered for an equal number of time steps until the
    % end of galileoRV
    if length(satT) > length(galileoRV)
        plot((satT(1:length(galileoRV))-satT(1))./24./3600,...
            (vecnorm(satRV(1:length(galileoRV),1:3)')-...
            vecnorm(galileoRV(1:3,:)))./au)
    
    % Otherwise, the opposite it true and galileoRV will be modeled until
    % the end of the items in satT
    else
        plot((satT-satT(1))./24./3600,...
            (vecnorm(satRV(:,1:3)')-...
            vecnorm(galileoRV(1:3,1:length(satRV))))./au)
    end
    xlabel("Time Since Departure (days)")
    ylabel("Distance (AU)")
    title("Difference in Position Magnitudes")
    
    % Plot of geneitc algorithm's trajectory versus Galileo, difference in 
    % position magnitude scaled to percent error
    figure
    if length(satT) > length(galileoRV)
        plot((satT(1:length(galileoRV))-satT(1))./24./3600,...
            (vecnorm(satRV(1:length(galileoRV),1:3)')-...
            vecnorm(galileoRV(1:3,:)))./galileoMag*100)
    else
        plot((satT-satT(1))./24./3600,...
            (vecnorm(satRV(:,1:3)')-...
            vecnorm(galileoRV(1:3,1:length(satRV))))./...
            galileoMag(1:length(satRV))*100)
    end
    xlabel("Time Since Departure (days)")
    ylabel("Percent Difference (%)")
    title("Percent Difference in GA Trajectory vs Galileo")
    grid on
end

cspice_kclear;
end

%% Subfunctions
function [deltaV penalty] = gadV(vIn, vOut, planetV, planetNum, tol)
% In this function, calculate the needed Delta-V maneuver to patch the
% inbound and outbound velocities together
%
% Inputs**
% vIn = inbound, heliocentric velocity vector
% vOut = outbound, heliocentric velocity vector
% planetV = gravity assist planet's heliocentric velocity vector
% planetNum = The planet number involed in the gravity assist
% tol = Function tolerance for determining the outbount eccentricities
%
% Outputs**
% deltaV = scalar value of the Delta-V needed to patch trajectories
%          together
% penalty = scalar sum of the penalities accured for bad trajectories 


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

% Compute the realtive inbound and outbound velocity vectors
vIn_Inf = vIn - planetV;
vOut_Inf = vOut - planetV;

% Compute the inbound and outbound semi-major axes
aIn = -mu(planetNum)/norm(vIn_Inf)^2;
aOut = -mu(planetNum)/norm(vOut_Inf)^2;

% solve for the outbound trajectory's eccentricity
eOutNew = eOutCalc(planetNum,vIn_Inf,vOut_Inf,tol);

% If the final solution is not a real number, then return a large Delta-V
% value and peanlty
if ~isreal(eOutNew)
    deltaV = 50;
    penalty = 50;
    return
end
    
% Compute the patched orbit periapsis radius
rP = aOut*(1 - eOutNew);

% Planetary Mean Radii
rMercury = 2439; % km
rVenus = 6051.8; % km
rEarth = 6371; % km
rMars = 3389.5; % km
rJupiter = 69911; % km
rSaturn = 58232; % km
rUranus = 25362; % km
rNeptune = 24622; % km
rPluto = 1188.3; % km
rPlanet = [rMercury, rVenus, rEarth, rMars, rJupiter, ...
           rSaturn, rUranus, rNeptune, rPluto];

% Compute the difference in velocities at the periapsis radius of the
% gravity assist
deltaV = abs(sqrt(norm(vIn_Inf)^2 + 2*mu(planetNum)/rP) - ...
            sqrt(norm(vOut_Inf)^2 + 2*mu(planetNum)/rP));

% List of planetary mass parameters
au = 1.495978e8; % km, 1 au
mercuryAU = 0.3871; % au
venusAU = 0.7233; % au
earthAU = 1; % au
marsAU = 1.5237; % au
jupiterAU = 5.2028; % au
saturnAU = 9.5388; % au
uranusAU = 19.1914; % au
neptuneAU = 30.0611; % au
plutoAU = 39.5294; % au
planetDistances = [mercuryAU, venusAU, earthAU, marsAU, jupiterAU, ...
saturnAU, uranusAU, neptuneAU, plutoAU];
% For each planet, compute the radius of the sphere of influence
for l = 1:length(mu)-1
    rSoi_planets(l) = (mu(l)/mu(10))^(2/5)*planetDistances(l)*au;
end

% For penalties, the first to check is if the flyby periapsis is below the
% safety threshold (1.05*rPlanet)
flybyCoeff = 1.05;
if rP < flybyCoeff*rPlanet(planetNum)
    % if the periapsis is bleow the threshold, assign a alrge penalty
    periapsisPenalty = -10*log10(rP/(flybyCoeff*rPlanet(planetNum)));
else
    % Otherwise do not add a penalty
    periapsisPenalty = 0;
end

% The second penalty to check for is a low energy transfer
E = (norm(0.9.*vIn)^2)/2 - mu(planetNum)/rSoi_planets(planetNum);
if E >= 0
    % If the energy is greater than or equal to 0 then it means the inbound
    % orbit is hyperbolic or parabolic
    lowEnergyPenalty = 0;
else
    % However, if the value is negative then it means that the inbound
    % orbit is elliptical.
    lowEnergyPenalty = 1/norm(vIn);
end

% Sum the penalties together for this gravity assist
penalty = periapsisPenalty + lowEnergyPenalty;
end

function eOutNew = eOutCalc(planetNum,vIn_Inf,vOut_Inf,tol)
% This subroutine is the newton-raphson iterative root finder to solve for
% the outbound eccentricity that helps patch the inbound and outbound
% trajectories

% Planetary Mass Parameters
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

% Get the inbound and outbound semi-major axes
aIn = -mu(planetNum)/norm(vIn_Inf)^2;
aOut = -mu(planetNum)/norm(vOut_Inf)^2;

% Compute the hyperbolic turning angle
turnAngle = acos(dot(vIn_Inf,vOut_Inf)/norm(vIn_Inf)/norm(vOut_Inf));

% Local function as a function of eOut. The function is set equal to 0
f = @(eOut) (aOut/aIn*(eOut-1))*sin(turnAngle - asin(1/eOut)) - 1;

% Take the derivative of the function with respect to eOut 
dfdeOut = @(eOut) (aOut/aIn*eOut - aOut/aIn + 1) ...
    * cos(turnAngle - asin(1/eOut))/(eOut^2*sqrt(1 - 1/eOut^2)) ...
    + aOut/aIn*sin(turnAngle - asin(1/eOut));


% Iteration counter
loop = 1;
% Starting guess for eOut
eOutOld = 1.5;
% Initial solution for Newton-Raphson
eOutNew = eOutOld - f(eOutOld)/dfdeOut(eOutOld);

% If the value ever drops below 1, this logical flips to try solving for
% eOut with a higher intial guess
low = false;

% While the difference between the new and old eOuts, perform the following 
while abs((eOutOld - eOutNew)/eOutOld) > tol
    % Update the values of eOld and eNew
    eOutOld = eOutNew;
    eOutNew = eOutOld - f(eOutOld)/dfdeOut(eOutOld);

    % If eNew is not real, that means eOld was less than 1.
    if ~isreal(eOutNew)
        if ~low
            % Try the root finder again but with a greater intial guess
            eOutNew = 10;
            % change value of low so that this will end if another nonreal
            % answer is computed
            low = true;
        elseif low
            % If eNew is imaginary for the second time, then just return
            % to the function call with the imaginary solution
            return
        end
        % reset the loop counter to give new root finder a change to find
        % the solution
        loop = 0;
    end
    
    % If the number of loops has reached 200, the value of eNew is less
    % than 0, or the value of eNew is NAN, perform the following
    if loop == 200 || eOutNew < 0 || isnan(eOutNew)
        % rather than continue mesing around with Newton-Raphson, use
        % MATLAB built-in non-linear function solver, fsolve
        options = optimset('Display','off');
        eOutNew = fsolve(f,1.5,options);
        
        % After solving for eNew, return to the function call
        return
    end

    % increase the loop counter
    loop = loop + 1;
end

end