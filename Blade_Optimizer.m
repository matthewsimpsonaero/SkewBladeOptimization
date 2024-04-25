% Matthew Simpson
% Turbine Blade Optimizer in Skew
% 1-27-2024
clc;clear, close all

% This code selects airfoils for the turbine optimization based on input
% paramters

addpath(genpath('./Initalization Functions'))
addpath(genpath('./Genetic Algorithm Functions'))
generate_plots = true; % Toggle this if you would like to show plots while the code is running

%% Input parameters

% These parameters are what are used when generating the Reynolds number
% and mach number of each blade element

V_inf = 5; % Inflow velocity (m/s)
Blade_length = 0.15; % (m)
skew_angle = 20; %Angle of flow relative to the turbine (deg)
blade_number = 3; % How many blades will the turbine have
TSR = (4*pi/blade_number)*1.25;  % Optimal TSR of the turbine from https://users.wpi.edu/~cfurlong/me3320/DProject/Ragheb_OptTipSpeedRatio2014.pdf
n = 10; %number of blade element segments
hub_radius = 0.0127; % How large the hub is. This will dictate where the first element starts.

%% Flow Paramters For density determination

% This can be changed as needed to model in other fluids such as water,
% right now it is set up for room temperature air at sea level

air_temp = 20; % celcius
atmos_pressure = 101325; % Sea level (Pa)
dynamic_viscosity_air = 1.8205e-5; % (kg/m*s)
gamma_air = 1.4; % ratio of specifc heats for air
R = 8.314459; % Gas constant (j/mol*K)
molar_mass_air = 28.9628; % (g/mol) https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html

Density_air = ComputeDensityAir(air_temp,atmos_pressure,R,molar_mass_air); % located in "Initalization functions"

%% Determination of velocity and direction

Element_location = linspace(0,Blade_length,n+1); % create equally spaced elements from hub to blade length, n+1 is because
% of how Qblade treats the first element as starting at 0.

omega = (TSR*V_inf)/ (Blade_length+hub_radius); %rotational rate of the turbine from TSR equation rearranged(rad/s)

Position = zeros(n+1,360); %preallocate to save time
Skew_Velocity = zeros(n+1,360);
for i = 1:n+1 % compute for each blade element
    for theta = 1:360 % compute for each turbine azimuth angle as it rotates
        Position(i,theta) = (Element_location(i)+hub_radius)*sind(skew_angle)*cosd(theta); % see AIAA Paper equation 5
        Skew_Velocity(i,theta) = -omega*(Element_location(i)+hub_radius)*sind(skew_angle)*sind(theta); % See AIAA Paper equation 6
    end
end

V_total_static = Skew_Velocity+(V_inf*cosd(skew_angle)); % total velocity without effects of the rotating turbine

if generate_plots
    GeneratePositionPlot(Position,n) % generates the position vs azimuth angle plot for the skew
    GenerateSkewVelocityPlot(Skew_Velocity,n) % generates the velocity vs azimtuh angle for the skew
    GenerateTotalAxialVelocityPlot(V_total_static,n) % generates the total velocity without effects of the rotating turbine vs azimuth angle
end

V_Y = zeros(n+1,360); %preallocate to save time
V_X = zeros(n+1,360);
Velocity_Direction = zeros(n+1,360);
theta = 1:360; % This needs to be here as the whole azimuth sweep is
% computed at once, eleminating the need for nested FOR loops
for i = 1:n+1
    % see slide 7 of confrence presentation for axis definition
    V_Y(i,1:360) = V_inf*cosd(skew_angle)-omega*(Element_location(i)+hub_radius)*sind(skew_angle)*sind(theta)*cosd(skew_angle); % see AIAA Paper equation 3 and 4
    V_X(i,1:360) = V_inf*sind(skew_angle)*sind(theta) + omega*(Element_location(i)+hub_radius);  % see AIAA Paper equation 2 and 6

    Velocity_Direction(i,1:360) = atan2d( V_Y(i,1:360), V_X(i,1:360)); % see slide 13 of confrence presentation
end

% Compute the average velocity direction. This is needed to optimize the
% turbine for a skew conditon, rather than for a single turbine azimuth
% angle. If we did not do this average, we could optimize the turbine for
% each angle from 0-360.
Average_Velocity_Direction = mean(Velocity_Direction,2)';

% plot the total velocity magnitude and direction vs. azimuth angle
if generate_plots
    PlotTotalVelocityMagnitude(V_Y,V_X,n)
    PlotTotalVelocityDirection(Velocity_Direction,n)
end

%% Compute the optimal chord distrobution per Betz equation

GuessCL = 1.5; % Provide an inital guess of the CL, this is done as to not
% have to solve for the chord distrobution iterativly.
roverR = (linspace(hub_radius/(Blade_length+hub_radius),1,n+1)); % compute
% the nondimensional element location by taking the element location
% divided by the total blade length

% compute the tip speed ratio at each location. This equation is long since
% the tip speed ratio does not start at 0, since there is a hub and the
% element closest to the hub has a non-zero rotational rate
LamdaR = hub_radius*TSR/((Blade_length+hub_radius)) + TSR/((Blade_length+hub_radius))*(roverR*Blade_length);

% uncomment lines 101 and 102 to remove the constant chord, and solve for
% the optimal chord using the Betz equations
coverR_wake = zeros(1,n+1); %preallocate to save time
coverR_no_wake = zeros(1,n+1);
for i = 1:n+1
    coverR_wake(i) = (1/3); % this 1/3 is to give the chord a constant length of 5 cm, compared to the length of 15 cm
    coverR_no_wake(i) = (1/3);
    % This equation is for when the wake is assumpted to be rotating. This
    % would be a more realistic model of the system
    %coverR_wake(i) = ((8*pi*roverR(i))/(blade_number*GuessCL))*(1-cosd(Average_Velocity_Direction(i))); %Wiley "Wind Energy Explained" equation 3-106
    % This equation is for when wake is modeled to be static. This would be
    % a less realistic model of the system
    %coverR_no_wake(i) = ((8*pi*roverR(i))*sind(Average_Velocity_Direction(i)))/(3*blade_number*GuessCL*LamdaR(i));% Wiley "Wind Energy Explained" 3-79
end

if generate_plots % generate a plot of the chord distrobution if needed
    GenerateChordPlot(roverR,coverR_wake,coverR_no_wake)
end

%% Compute the average velocity and average mach number at each element

Average_Velocity = zeros(1,n+1); %preallocate to save time
Average_Mach = zeros(1,n+1);
for i = 1:n+1
    % Find the average velocity by finding the mean of the vector magnitude
    % of V_X and V_Y
    Average_Velocity(i) = mean(sqrt(V_X(i,:).^2+V_Y(i,:).^2));
    % Find the average mach by taking the vector magnitude of V_X and V_Y
    % and dividing by the speed of sound formula at room air temperature
    Average_Mach(i) = mean(sqrt(V_X(i,:).^2+V_Y(i,:).^2) ./ sqrt(gamma_air*(R*1000/molar_mass_air)*(air_temp+273.15)));
end

%% Compute thee Reynolds Number

Reynolds_Number = zeros(n+1,360); %preallocate to save time
for i = 1:n+1
    % Compute the reynolds number at each azimuth angle using the reynolds 
    % number formula
    Reynolds_Number(i,1:360) = (Density_air*sqrt(V_X(i,:).^2+V_Y(i,:).^2)*(coverR_wake(i)*Blade_length))/dynamic_viscosity_air;
end

Average_Reynolds_Number = mean(Reynolds_Number,2)'; % Find the average 
% reynolds number

if generate_plots % plot reynolds number vs turbine azimuth angle
    PlotReynoldsNumber(Reynolds_Number,n)
end

%% Xfoil Selection Using a Genetic Algorithm (Setup)

fid2 = fopen( 'results.txt', 'wt' ); % results will store the final GA output once all the generations have run
fid3 = fopen( 'Generational_Best.txt', 'wt' ); % generational best will store the GA outputs for each generation 
fprintf(fid2, 'Element Number\t\tAirfoil\t\tint(Cl/Cd)\t\tReynolds#\n'); % print headers
fprintf(fid3, 'Element Number\t\tAirfoil\t\tint(Cl/Cd)\t\tReynolds#\n');

% This algorithm selects from NACA 4 series airfoils. There are 9*9*40 =
% 3240 possible cobinations of airfoils

starting_population = 600; % select how many airfoils will be generated as the inital population
num_generations = 5; % select the number of generations
num_parents = [200,100,50,20]; % select how many parents will be selected. 
% These parents will have the highest integral(CL/CD) across the angle of attack sweep
num_offspring = [200,100,50,20]; % select how many offspring will be generated
carryover = 5; % take the top 5 parents and move them to the next generation. This ensures that best airfoils are not left behind

course_xfoil = 40; % Number of angle of attack points that will be ran 
points_course = linspace(-20,20,course_xfoil); % The angle of attck will be from -20 deg to 20 deg. This can be changed

%% Xfoil Selection Using a Genetic Algorithm (Run)

% By the end of the run, generational_best.txt and results.txt will be
% fully populated. During the run, an actions.txt file will be generated.
% This is telling Xfoil what to do. There will also be a tp......txt file
% generated, this stores the airfoil points to feed to xfoil

for element_num = 1:n+1
    % select the reynolds number and mach number for the specific element
    % that is being run
    Element_Reynolds = Average_Reynolds_Number(element_num);
    Element_Mach = Average_Mach(element_num);
    
    % generate the inital population. See the function for details
    Population = generate_inital_population(starting_population);
    
    for gen = 1:num_generations % iterate through each generation
        for i = 1:length(Population) % iterate through each airfoil in the population
            NACA_Number = Population{i}; % get the airfoil from the population
            try % sometimes xfoil fails, using a try/catch allows the code to continue
                
                fprintf('G%d:(%d/%d)Running NACA%s\n',gen,i,length(Population),NACA_Number) % print status to command window

                [Polar] = Airfoil_Runner(NACA_Number,Element_Reynolds,Element_Mach,points_course); %https://www.mathworks.com/matlabcentral/fileexchange/30478-rafael-aero-xfoilinterface
                LiftoverDragint(i) = trapz(Polar.Alpha,(Polar.CL./Polar.CD)); % compute the integral(CL/CD) across the angle of attack sweep 
            catch
                LiftoverDragint(i) = 0; % if xfoil fails, save as 0. This will allow the airfoil to be discarded on the next generation
            end
        end
        if gen~=num_generations % if the final generation is not reached
        [sorted_values, sorted_indexes] = sort(LiftoverDragint, 'descend'); % sort the integral(CL/CD) by the output value
        top_values = sorted_values(1:num_parents(gen)); % select the top values from 1->num parents. This will generate parent airfoils. This is not used but it is stored for metrics
        indexes = sorted_indexes(1:num_parents(gen)); % store the indexes of the parents to feed into the parents cell array
        Parents = {};
        for i = 1:num_parents(gen)
            Parents{i} = Population(indexes(i)); % pull the NACA 4-series number of the parents from the total population
        end
        fprintf(fid3, '%d\t\t\t\t\tNACA%s\t%.3f\t\t\t%d\n',element_num,Population{sorted_indexes(1)},sorted_values(1),Element_Reynolds); % print the generational best to the file
        Population = generate_offspring(Parents,num_offspring(gen)); % see function for details
        Population = [Population , [Parents{1:carryover}]];  % carry over best performing parents into the next generation
        clear LiftoverDragint % clear the metric variable for use in the next generation. As the next generation will be smaller, indexes greater than the number of offspirng need to be removed
        
        else
        [sorted_values, sorted_indexes] = sort(LiftoverDragint, 'descend'); % sort the integral(CL/CD) by the output value
        indexes = sorted_indexes(1); % get the top airfoil from the sort
        end
    end
    
    % run a fine Xfoil Analysis of the best airfoil to get a higher fidelity AOA value
    [Maxlift,bestAOA,Lift_curve,Drag_curve,AOA_curve,maxLoverDint] = RunFineXfoil(Population{indexes(1)},Element_Reynolds,Element_Mach);
    
    fprintf('Best Airfoil: %s\n',Population{indexes(1)}) % print status to command window
    fprintf(fid2, '%d\t\t\t\t\tNACA%s\t%.3f\t\t\t%d\n',element_num,Population{indexes(1)},maxLoverDint,Element_Reynolds); % print the generational best airfoil to the file using the output values from runfilexfoil function
end


