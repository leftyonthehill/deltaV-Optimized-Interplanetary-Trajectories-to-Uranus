%% Multiple Gravity Assist Trajectory Design
% This routine uses a genetic algorithm to find optimized solutions that
% minimize the total required Delta-V while using multiple planetary
% gravity assists. Refer to mgaDeisgnV2.m, mgaTOFDesignV2.m, or
% mgaPlotterV2.m for for information.
%
% For reference, any use of the planets will be referred to in their order
% from the sun. Venus will be represented by a "2" and Jupiter will be
% represented by a "5", for example.
%
% **WARNING** The default option for the genetic algorithm function in this
% script is to use parallel processing to compute the trajectories. If you
% do not want to use this setting, please update line 259.
%
% **WARNING** With a large number of output trajectories, the computation
% time will take many hours and large portions of RAM. From testing, this
% script when run in parallel to compute 60 trajectories took over 12 hours
% and more than 15 GM of RAM

% User variables
mission = "Galileo"; % Select either "Galileo" "Uranus"
numTrajectories = 1;
toPlot = 1;
% Galileo Search Conditions
if lower(mission) == lower("Galileo")
    % For the Galileo recreation, only one trajectory will be found. The
    % cose however may need to be run multiple times to get a valid
    % recreation
    numTrajectories = 1;
    toPlot = 1;

    % Starting plant
    planetStart = 3;
    % Parking orbit radius
    rDepart = 6378.14+296;

    % Target planet
    planetTarget = 5;
    % Target orbit shape
    targetOrbitShape = [286000 0.998];
    
    % Inertial Upper Stage injection energy
    departC3 = 17;
    
    % Launch window start date and time
    startDate = "Oct 01, 1989";
    startTime = "00:00:00.000";
    
    % Converts input time to a time understandable by the cspice toolkit
    t0 = datetime(startDate+" "+startTime, ...
        'InputFormat',"MMM dd, yyyy HH:mm:ss.SSS");
    t0 = datetime(t0, ...
        'Format',"MMM dd, yyyy HH:mm:ss.SSS");
    
    % Dates throughout this script, mgaDeisgnV2.m, mgaTOFDesignV2.m, and
    % mgaPlotterV2.m use Julian days to assist in the randomization of the
    % genetic algorithm decision variables
    t0Departure = juliandate(t0);
    
    % End of the launch window search
    endDate = "Nov 30, 1989";
    endTime = "00:00:00.000";
    
    % Converts input time to a time understandable by the cspice toolkit
    tf = datetime(endDate+" "+ endTime, ...
        'InputFormat',"MMM dd, yyyy HH:mm:ss.SSS");
    tf = datetime(tf, ...
        'Format',"MMM dd, yyyy HH:mm:ss.SSS");
    
    tfDeparture = juliandate(tf);
    
    % Arrays to store the upper and lower bound constraints for the genetic
    % algorithm decision variables. The order of the data in these arrays
    % must be as follows:
    % 1) Launch day (in Julian days)
    % 2) Number of gravity assists in the trajectory
    % 3) Order of gravity assist (will be more than the number of allowed 
    %    gravity assists)
    % 4) Times of flight between gravity assists
    % 5) Final trajectory leg time of flight
    % 6) B-plane angles for when Gooding's lambert's problem solver cannot
    %    find a solution
    minConditions = [];
    maxConditions = [];
    
    % Assigning the upper and lower bound of the launch window
    launchWindowSearchOpen = t0Departure;
    launchWindowSearchClose = tfDeparture;
    minConditions = [minConditions launchWindowSearchOpen];
    maxConditions = [maxConditions launchWindowSearchClose];
    
    % For the Galileo recreation, these values are fixed but this is the
    % upper and lower bound assignment for the number of gravity assists
    maxN = 3;
    minFlyBy = 3;
    minConditions = [minConditions minFlyBy];
    maxConditions = [maxConditions maxFlyBy];
    
    % A counting variable useful for later
    maxFlyBy = maxN;

    % Assignment of the upper and lower bounds of which planets can
    % contribute gravity assists
    minConditions = [minConditions 2 3 3];
    maxConditions = [maxConditions 2 3 3];
    
    % The commented out section was original test to recreat Galileo but
    % no good solution emerged. Later tightened the bounds closer to
    % Galileo's true flight times for a better solution.
    % gaLegMinLength = 100;
    % gaLegMaxLength = 1200;
    % for i = 1:maxN
    %     minConditions = [minConditions gaLegMinLength];
    %     maxConditions = [maxConditions gaLegMaxLength];
    % end
    % 
    % gaFinalLegMinLength = 100;
    % gaFinalLegMaxLength = 1200;
    % minConditions = [minConditions gaFinalLegMinLength];
    % maxConditions = [maxConditions gaFinalLegMaxLength];
    
    % From L. A. D’Amario, L. E. Bright, and A. A. Wolf, “Galileo 
    %      trajectory design: The Galileo mission,” Space science 
    %      reviews, vol. 60, no. 1–4, pp. 23–78, 1992.
    %
    % Assignment of the upper and lower bounds for the times of flight
    % Galileo Key Dates
    % Departure: 10/18/1989
    % Venus Close approach: 02/10/1990 (115)
    % Earth 1 Close approach: 12/08/1990 (301)
    % Earth 2 Close approach: 12/08/1992 (731)
    % Jupiter Arrival (closest approach): 12/07/1995 (1094)
    minConditions = [minConditions 105 291 721 1064];
    maxConditions = [maxConditions 125 311 741 1124];


% Uranus Search Conditions
elseif lower(mission) == lower("Uranus")
    % Starting planet
    planetStart = 3;
    % Starting parking orbit
    rDepart = 6378.14+250;
    
    % Target planet
    planetTarget = 7;
    % Target orbit shape
    targetOrbitShape = [300000 0.998];
    
    % Inertial Upper Stage injection energy
    departC3 = 18;
    
    % Start of the search window
    startDate = "Jan 01, 2030";
    startTime = "00:00:00.000";
    
    % Converts input time to a time understandable by the cspice toolkit
    t0 = datetime(startDate+" "+startTime, ...
        'InputFormat',"MMM dd, yyyy HH:mm:ss.SSS");
    t0 = datetime(t0, ...
        'Format',"MMM dd, yyyy HH:mm:ss.SSS");
    
    t0Departure = juliandate(t0);
    
    % End of the search window
    endDate = "Dec 31, 2034";
    endTime = "00:00:00.000";
    
    % Converts input time to a time understandable by the cspice toolkit
    tf = datetime(endDate+" "+ endTime, ...
        'InputFormat',"MMM dd, yyyy HH:mm:ss.SSS");
    tf = datetime(tf, ...
        'Format',"MMM dd, yyyy HH:mm:ss.SSS");
    
    tfDeparture = juliandate(tf);
    
    % Arrays to store the upper and lower bound constraints for the genetic
    % algorithm decision variables. The order of the data in these arrays
    % must be as follows:
    % 1) Launch day (in Julian days)
    % 2) Number of gravity assists in the trajectory
    % 3) Order of gravity assist (will be more than the number of allowed 
    %    gravity assists)
    % 4) Times of flight between gravity assists
    % 5) Final trajectory leg time of flight
    % 6) B-plane angles for when Gooding's lambert's problem solver cannot
    %    find a solution
    minConditions = [];
    maxConditions = [];
    
    % Assignment of the upper and lower bounds of the search window
    launchWindowSearchOpen = t0Departure;
    launchWindowSearchClose = tfDeparture;
    minConditions = [minConditions launchWindowSearchOpen];
    maxConditions = [maxConditions launchWindowSearchClose];
    
    % The upper and lower bound on the number of gravity assists
    minFlyBy = 2;
    maxFlyBy = 10;
    minConditions = [minConditions minFlyBy];
    maxConditions = [maxConditions maxFlyBy];
    
    % A counting variable useful for later
    maxN = maxFlyBy;
    
    % The upper and lower bound for each possible gravity assist
    minFlyByPlanet = 1;
    maxFlyByPlanet = 5;
    for i = 1:maxN
        minConditions = [minConditions minFlyByPlanet];
        maxConditions = [maxConditions maxFlyByPlanet];
    end
    
    % The upper and lower bounds on the times of flight for each of the
    % intermediate legs of the trajectory
    gaLegMinLength = 50;
    gaLegMaxLength = 2000;
    for i = 1:maxN
        minConditions = [minConditions gaLegMinLength];
        maxConditions = [maxConditions gaLegMaxLength];
    end
    
    % Upper and lower bound on the final leg's time of flight
    gaFinalLegMinLength = 100;
    gaFinalLegMaxLength = 6000;
    minConditions = [minConditions gaFinalLegMinLength];
    maxConditions = [maxConditions gaFinalLegMaxLength];

% If the mission presented is not Galileo or Uranus
else
    disp("Invalid mission selection. Please select either Galileo " + ...
        "or Uranus.\nOr, design custom mission parameters in the " + ...
        "else block starting on line 224.")
end

% Should Gooding's Lambert solver fail to produce a trajectory, these are
% b-plane angles used in the computation of a resonance orbit
minAngle = -2*pi;
maxAngle = 2*pi;
minConditions = [minConditions minAngle minAngle];
maxConditions = [maxConditions maxAngle maxAngle];


%% Computation of N Gravity Assists Using the Genetic Algorithm

% Start timer
tic;

% Storage of identified trajectories
optimalTrajectory = [];

% number of found trajectories
count = 1;

% While there are less computed trajectories than desired, continue to
% perform the following...
while count <= numTrajectories
    % Options for the genetic algorithm
    gaOption = optimoptions("ga", ...
            ... % Maximum generations without improvement
                MaxStallGenerations=25, ... 
            ... % Change in tolerance that defines improvement
                FunctionTolerance=1e-5, ...
            ... % Parallel computation status
                UseParallel=true, ...       
            ... % Output the iteration information
                Display="iter");
    
    % Genetic algorithm function call. Collect the optimize decision
    % variables and their cost function score
    [x,val] = ga( ...
                @(x)mgaDesignV2( ...      % The fitness/cost function
                    planetStart, ...      % Starting planet
                    rDepart, ...          % Parking orbit
                    departC3, ...         % IUS injection energy
                    planetTarget, ...     % Target planet
                    targetOrbitShape, ... % Target orbit shape
                    maxN, ...             % Maximum gravity assists
                    optimalTrajectory, ...% Matrix of identified orbits
                    x), ...               % Array to optimize
                length(minConditions), ...% Number of variables to optimize
                [],[],[],[], ...          
                minConditions, ...        % Array of lower bounds
                maxConditions, ...        % Array of upper bounds
                [], ...
                2:2+maxN,gaOption);       % Variables that are integers
    
    % If there are saved trajectories, perform the following...
    if ~isempty(optimalTrajectory)
        % Get the number of optimized flybys
        numFlyBys = x(2);
        % Get the order of the gravity assists
        flightOrder = x((1:numFlyBys) + 2);

        % If there exists trajectories with the same number of graivty
        % assists, get their indexes
        idx1 = find(optimalTrajectory(:,3) == numFlyBys);

        % In the list of idx1, if any listing has the same order of gravity
        % assists, get their index
        idx2 = find(all( ...
            optimalTrajectory(idx1,((1:numFlyBys) + 3)) == flightOrder,2));

        % If idx2 comes back empty, the latest calculated trajectory is not
        % a duplicate
        if isempty(idx2)
            % Save the cost function score, departure time, and all other
            % variables except the number of gravity assists and the 
            % gravity assist order
            testTrajectory = [val x(1) x((3+maxN):end)];
            
            % Now the goal is to improve the initial solution by running it
            % through the genetic algorithm additional times. However,
            % since the route is fixed, the number of variables has
            % changed. The size of the upper and lower bound array need to
            % change befire the genetic algoritm can be recalled

            % Save the trajectory lower bound times of flight
            tofs = minConditions((maxN+3):(2*maxN+2));
            % Get the terminal leg's lower bound time of flight
            tf = minConditions(2*maxN+3);
            % Save an array of all trajectory lower bound times
            testMinBound = [minConditions(1);
                            tofs;
                            tf;
                            minConditions(end-1:end)]';
    
            % Save the trajectory upper bound times of flight
            tofs = maxConditions((maxN+3):(2*maxN+2));
            % Save the terminal leg's upper bound time of flight
            tf = maxConditions(2*maxN+3);
            % Save an array of all trajectory upper bound times
            testMaxBound = [maxConditions(1);
                            tofs;
                            tf;
                            maxConditions(end-1:end)]';
    
            
            % Save the last cost value as the current best value
            valBest = val;
            % Loop counter
            innerloop = 1;
            % While less than or equal to 5 loops, perform the following
            while innerloop <= 5
                % This loop is designed to improve the initial solution
                % found with the genetic algorithm. The found trajecotry
                % will improve usually at least one but sometimes not at
                % all.
                
                % Secondary genetic algorithm function call. All inputs are
                % the same except that there are now no interger values to
                % consider. Also, note that the fitness function has
                % changed. This modification of the original fitness
                % function is designed to optimize only the times of flight
                [xTest,valTest] = ga( ...
                    @(x)mgaTOFDesignV2( ...
                        planetStart, ...
                        rDepart, ...
                        departC3, ...
                        planetTarget, ...
                        targetOrbitShape, ...
                        maxN,numFlyBys, ...
                        flightOrder, ...
                        testTrajectory, ...
                        x), ...
                    length(testMinBound), ...
                    [],[],[],[], ...
                    testMinBound, ...
                    testMaxBound, ...
                    [], ...
                    [],gaOption);

                % If the output of the secondary genetic algorithm is
                % better than the previous cost value, then override the
                % best cost with the new value
                if valTest < valBest
                    testTrajectory = [valTest xTest];
                    valBest = valTest;
                end

                % Increase the loop
                innerloop = innerloop + 1;
            end
    
            % After reoptimizing the current trajectory, save its results
            % in optimalTrajectory and print to the terminal that a new
            % trajectory was identified
            planets = zeros(1,maxN);
            planets(1:numFlyBys) = flightOrder;
            optimalTrajectory(count,:) = [testTrajectory(1:2);
                                          numFlyBys;
                                          planets;
                                          testTrajectory(3:end)]';


            str = "New fly-by order: " + planetStart + "-";
            for i = 1:numFlyBys
                str = str + flightOrder(i) + "-";
            end
            str = str + planetTarget;
            disp(str)
            count = count + 1;

        % Although incredibly unlikely, if the genetic algorithm converges
        % onto a traejctory that improves ones of the saved routes, then
        % override its location within optimalTrajectories
        elseif optimalTrajectory(idx2,1) > val
            optimalTrajectory(idx2,:) = [val x];
        % In a even less likely scenario, if the genetic algorithm
        % converges onto a solution that is worse that what is known, then
        % do nothing and continue the loop. Do not increase the counter as
        % a new trajectory was not identified
        else
            continue
        end
    
    % If there are currently no saved trajectories, then perform the
    % following...
    else
        % Get the number of optimized flybys
        numFlyBys = x(2);
        % Get the order of the gravity assists
        flightOrder = x((1:numFlyBys) + 2);

        % Save the cost function score, departure time, and all other
        % variables except the number of gravity assists and the 
        % gravity assist order
        testTrajectory = [val x(1) x((3+maxN):end)];
        
        % Now the goal is to improve the initial solution by running it
        % through the genetic algorithm additional times. However,
        % since the route is fixed, the number of variables has
        % changed. The size of the upper and lower bound array need to
        % change befire the genetic algoritm can be recalled

        % Save the trajectory lower bound times of flight
        tofs = minConditions((maxN+3):(2*maxN+2));
        % Get the terminal leg's lower bound time of flight
        tf = minConditions(2*maxN+3);
        % Save an array of all trajectory lower bound times
        testMinBound = [minConditions(1)...
                        tofs...
                        tf...
                        minConditions(end-1:end)];

        % Save the trajectory upper bound times of flight
        tofs = maxConditions((maxN+3):(2*maxN+2));
        % Save the terminal leg's upper bound time of flight
        tf = maxConditions(2*maxN+3);
        % Save an array of all trajectory upper bound times
        testMaxBound = [maxConditions(1)...
                        tofs...
                        tf...
                        maxConditions(end-1:end)]';
    
        % Save the last cost value as the current best value
        valBest = val;
        % Loop counter
        innerloop = 1;
        % While less than or equal to 5 loops, perform the following
        while innerloop <= 5
            % This loop is designed to improve the initial solution
            % found with the genetic algorithm. The found trajecotry
            % will improve usually at least one but sometimes not at
            % all.
            
            % Secondary genetic algorithm function call. All inputs are
            % the same except that there are now no interger values to
            % consider. Also, note that the fitness function has
            % changed. This modification of the original fitness
            % function is designed to optimize only the times of flight
            [xTest,valTest] = ga( ...
                @(x)mgaTOFDesignV2( ...
                    planetStart, ...
                    rDepart, ...
                    departC3, ...
                    planetTarget, ...
                    targetOrbitShape, ...
                    maxN, ...
                    numFlyBys, ...
                    flightOrder, ...
                    testTrajectory, ...
                    x), ...
                length(testMinBound), ...
                [],[],[],[], ...
                testMinBound, ...
                testMaxBound, ...
                [], ...
                [],gaOption);

            % If the output of the secondary genetic algorithm is
            % better than the previous cost value, then override the
            % best cost with the new value
            if valTest < valBest
                testTrajectory = [valTest xTest];
                valBest = valTest;
            end
            % Increase the counter
            innerloop = innerloop + 1;
        end

        % After reoptimizing the current trajectory, save its results
        % in optimalTrajectory and print to the terminal that a new
        % trajectory was identified
        planets = zeros(1,maxN);
        planets(1:numFlyBys) = flightOrder;
        optimalTrajectory(count,:) = [testTrajectory(1:2) ...
            numFlyBys planets testTrajectory(3:end)];

        
        str = "New fly-by order: " + planetStart + "-";
        for i = 1:numFlyBys
            str = str + flightOrder(i) + "-";
        end
        str = str + planetTarget;
        disp(str)
        count = count + 1;
    end
end

toc;
%% Data colelction for the report

% If the mission type if Uranus, get the key trajectory dates and gravity
% assist order for the top 13 trajectories. Although the list in the report
% details 10, the 7th, 8th, and 12th best trajectories directly intercept
% Venus are not valid trajectories
if lower(mission) == lower("Uranus") && numTrajectories > 1
    temp = sortrows(optimalTrajectory,1);
else
    temp = [val x];
end

for i = 1:toPlot
    % Get the number of flybys
    numFlyBys = temp(i,3);
    % Get the gravity assist order
    flyByOrder = temp(i,4:3+numFlyBys);

    % Small loop to print out the order of the trajectory
    str = "E-";
    for j = 1:length(flyByOrder)
        if flyByOrder(j) == 1
            str = str + "Me-";
        elseif flyByOrder(j) == 2
            str = str + "V-";
        elseif flyByOrder(j) == 3
            str = str + "E-";
        elseif flyByOrder(j) == 4
            str = str + "M-";
        elseif flyByOrder(j) == 5
            str = str + "J-";
        end
    end

    str = str + "U";
    disp(str)

    % Print the cost of the trajectory
    temp(i,1)

    % Print the launch date
    date1 = datetime(temp(i,2),'ConvertFrom','juliandate','Format', ...
        "MMM dd, yyyy HH:mm:ss.SSS")
    numFlyBys = temp(i,3);
    
    % Print the dates of closest approach for each gravity assist
    for j = 1:length(flyByOrder)
        date = datetime(temp(i,2)+sum(temp(i,maxN+4:maxN+3+j)),...
            'ConvertFrom','juliandate','Format', ...
            "MMM dd, yyyy HH:mm:ss.SSS")
    end

    % Print the arrival date
    date2 = datetime(...
        temp(i,2)+sum(temp(i,maxN+4:maxN+3+numFlyBys))+temp(i,2*maxN+4),...
        'ConvertFrom','juliandate','Format', ...
        "MMM dd, yyyy HH:mm:ss.SSS")

    % Print the time span from launch to arrival
    span = years(date2-date1)
    
    % Using a modified version of the genetic algorithm fitness
    % function, plot the trajectories
    C = mgaPlotterV2( ...
        mission, ...
        planetStart, ...
        rDepart, ...
        departC3, ...
        planetTarget, ...
        targetOrbitShape, ...
        maxN, ...
        temp(i,2:end));
end


