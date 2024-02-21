% Matthew Simpson
% Turbine Blade Optimizer in Skew
% 1-27-2024
clc;clear;close all
addpath(genpath('./Initalization Functions'))
addpath(genpath('./Genetic Algorithm Functions'))
generate_plots = true;

%% Input parameters

V_inf = 5; %m/s
Blade_length = 0.15; %m
skew_angle = 30; %deg
blade_number = 3;
TSR = (4*pi/blade_number)*1.25;  %https://users.wpi.edu/~cfurlong/me3320/DProject/Ragheb_OptTipSpeedRatio2014.pdf
n = 10; %number of blade element segments
hub_radius = 0.0127; %meters

%% Flow Paramters

air_temp = 20; %celcius
atmos_pressure = 101325; %Pa
dynamic_viscosity_air = 1.8205e-5; % kg/m*s
gamma_air = 1.4;
R = 8.314459; %j/mol*K
molar_mass_air = 28.9628; %g/mol

Density_air = ComputeDensityAir(air_temp,atmos_pressure,R,molar_mass_air);

%% Determination of velocity and direction

Element_location = linspace(Blade_length/n,Blade_length,n);
omega = (TSR*V_inf)/ (Blade_length+hub_radius);

for i = 1:n
    for theta = 1:360
        Position(i,theta) = (Element_location(i)+hub_radius)*sind(skew_angle)*cosd(theta);
        Skew_Velocity(i,theta) = -omega*(Element_location(i)+hub_radius)*sind(skew_angle)*sind(theta);
    end
end

V_total_static = Skew_Velocity+(V_inf*cosd(skew_angle));

if generate_plots
    GeneratePositionPlot(Position,n)
    GenerateSkewVelocityPlot(Skew_Velocity,n)
    GenerateTotalAxialVelocityPlot(V_total_static,n)
end

theta = 1:360;
for i = 1:n
    V_perpindicular(i,1:360) = V_inf*(cosd(theta)-sind(theta)*cosd(skew_angle));
    V_parallel(i,1:360) = V_inf*sind(skew_angle)*sind(theta) + omega*(Element_location(i)+hub_radius);
    Velocity_Direction(i,1:360) = atan2d( V_perpindicular(i,1:360), V_parallel(i,1:360));
end

Average_Velocity_Direction = mean(Velocity_Direction');

if generate_plots
fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,sqrt(V_parallel(i,:).^2+V_perpindicular(i,:) .^2))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Seen Velocity (m/s)')
title('Seen Velocity vs. Rotation Angle for Skewed Turbine')

fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,Velocity_Direction(i,:))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Velocity Direction (°)')
title('Velocity Direction vs. Rotation Angle for Skewed Turbine')
end

% Compute optimal chord with wake rotation Guessing CL
GuessCL = 1.5;
roverR = (linspace(1/n,1,n));
LamdaR = TSR*roverR;

for i = 1:n
    coverR_wake(i) = ((8*pi*roverR(i))/(blade_number*GuessCL))*(1-cosd(Average_Velocity_Direction(i))); %Wiley equation 3-106
    coverR_no_wake(i) = (8*pi*(Element_location(i)+hub_radius)*sind(Average_Velocity_Direction(i)))/(3*blade_number*GuessCL*LamdaR(i));% Wiley 3-79
end

if generate_plots
GenerateChordPlot(roverR,coverR_wake)
end
% compute solidarity
R_values = Element_location+hub_radius;
R_total = Blade_length+hub_radius;
integral_C = trapz(R_values, coverR_wake*(R_total));
sigma = ((integral_C*blade_number) +(pi*hub_radius^2)) / (pi * R_total^2); %Wiley 3.107

% Compute axial and tangential induction factors

for i = 1:n
aprime(i) = 1/(((4*cosd(Average_Velocity_Direction(i)))/(sigma*GuessCL))-1); % Wiley 3.89
a(i) = -((aprime(i)/((sigma*GuessCL*cosd(Average_Velocity_Direction(i)))/(4*LamdaR(i)*sind(Average_Velocity_Direction(i)))))-1); %wiley 3-83
end

% Compute Velocity using the induction factors

for i = 1:n
V_parallel_corrected(i,:) = (V_parallel(i,:).*(1+aprime(i)));
V_perpindicular_corrected(i,:) = (V_perpindicular(i,:).*(1-a(i)));
Velocity_Direction_corrected(i,1:360) = atan2d(V_perpindicular_corrected(i,:), V_parallel_corrected(i,:));
end

if generate_plots
fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,sqrt(V_parallel_corrected(i,:).^2+V_perpindicular_corrected(i,:) .^2))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Seen Velocity (m/s)')
title('Seen Velocity vs. Rotation Angle for Skewed Turbine')

fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,Velocity_Direction_corrected(i,:))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Velocity Direction (°)')
title('Velocity Direction vs. Rotation Angle for Skewed Turbine')
end

for i = 1:n
Direction_range(i) = max(Velocity_Direction_corrected(i,:))-min(Velocity_Direction_corrected(i,:));
Direction_min(i) = min(Velocity_Direction_corrected(i,:));
end


for i = 1:n
Average_Velocity(i) = mean(sqrt(V_parallel_corrected(i,:).^2+V_perpindicular_corrected(i,:).^2));
end

for i = 1:n
Average_Mach(i) = mean(sqrt(V_parallel_corrected(i,:).^2+V_perpindicular_corrected(i,:).^2) ./ sqrt(gamma_air*(R*1000/molar_mass_air)*(air_temp+273.15)));
end

%% Calculation of reynolds number

for i = 1:n
Reynolds_Number(i,1:360) = (Density_air*sqrt(V_parallel_corrected(i,:).^2+V_perpindicular_corrected(i,:).^2)*(coverR_wake(i)*Blade_length))/dynamic_viscosity_air;
end

if generate_plots
fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,Reynolds_Number(i,:))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Reynolds Number')
title('Reynolds Number vs. Rotation Angle for Skewed Turbine')
end

Average_Reynolds_Number = mean(Reynolds_Number');


%% Xfoil Selection Using a Genetic Algorithm

% fid2 = fopen( 'results.txt', 'wt' );
% fid3 = fopen( 'Generational_Best.txt', 'wt' );
% fprintf(fid2, 'Element Number\t\tAirfoil\t\tintCLoverCD\t\tAOA\n');
% fprintf(fid3, 'Element Number\t\tAirfoil\t\tint(l/d)\t\tReynolds#\n');
% 
% starting_population = 1200; % 50 percent of the total subset
% num_generations = 2;
% num_parents = starting_population/10;
% num_offspring = starting_population/10;
% 
% course_xfoil = 40; % number of points
% points_course = linspace(-20,20,course_xfoil);
% for element_num = 1:n
% 
%     Element_Reynolds = Average_Reynolds_Number(element_num);
%     Element_Mach = Average_Mach(element_num);
% 
%     Population = generate_inital_population(starting_population);
%     for gen = 1:num_generations
%         for i = 1:length(Population)
%             NACA_Number = Population{i};
%             try
%                 fprintf('G%d:(%d/%d)Running NACA%s\n',gen,i,length(Population),NACA_Number)
% 
%                  %[Polar] = Airfoil_Runner(NACA_Number,1,Element_Mach,points_course);
%                 [Polar] = Airfoil_Runner(NACA_Number,Element_Reynolds,Element_Mach,points_course);
%                 LiftoverDragint(i) = trapz(Polar.Alpha,(Polar.CL./Polar.CD));
%             catch
%                 LiftoverDragint(i) = 0;
%             end
%         end
%         [sorted_values, sorted_indexes] = sort(LiftoverDragint, 'descend');
%         top_values = sorted_values(1:num_parents);
%         indexes = sorted_indexes(1:num_parents);
%         Parents = {};
%         for i = 1:num_parents
%             Parents{i} = Population(indexes(i));
%         end
%         fprintf(fid3, '%d\t\tNACA%s\t\t%.3f\t\t%d\n',element_num,Population{sorted_indexes(1)},sorted_values(1),Element_Reynolds);  
%         if gen~=num_generations
%             Population = generate_offspring(Parents,num_offspring);
%             Population = [Population , [Parents{1:20}]];
%             clear LiftoverDragint
%         end
%     end
%     
%     % run a fine Xfoil Analysis to get a higher fidelity AOA value
%     [Maxlift,bestAOA,Lift_curve,Drag_curve,AOA_curve,maxLoverDint] = RunFineXfoil(Population{indexes(1)},Element_Reynolds,Element_Mach);
%     
%     fprintf('Best Airfoil: %s\n',Population{indexes(1)})
%     fprintf(fid2, '%d\t\tNACA%s\t\t%.3f\t\t%d\n',element_num,Population{indexes(1)},maxLoverDint,Element_Reynolds);
% end


%% Load best airfoils

fileID = fopen('Emelia Results/maximization/results.txt', 'r');
fgetl(fileID);
C = textscan(fileID, '%s', 'Delimiter', '\n');
dataCellArray = C{1};
fclose(fileID);

for i = 1:8
    temp = dataCellArray{i};
    if i<10
        NACAnum(i) = str2double(temp(8:11));
    else
        NACAnum(i) = str2double(temp(9:12));
    end
end

liftofElement = {};
AOAofElement ={};
Lift_curves = {};
Drag_curves = {};
AOA_curves = {};

for i = 1:8
[Maxlift,bestAOA,Lift_curve,Drag_curve,AOA_curve] = RunFineXfoil(num2str(NACAnum(i)),Average_Reynolds_Number(i),Average_Mach(i));

liftofElement{i} = Maxlift;
AOAofElement{i} = bestAOA;
Lift_curves{i} = Lift_curve;
Drag_curves{i} = Drag_curve;
AOA_curves{i} = AOA_curve;
end


%% Determine twist
Qblade_pitch = [48.1,26.4,16.61,11.3,6.098,3.917,1.353,.178,-0.736];

for jj = 1:8
[alpha_ext, CL_ext, CD_ext] = viterna_extrapolation(AOA_curves{jj}, Lift_curves{jj}, Drag_curves{jj}); %perform vinerna expansion 

% Define the known values (replace these placeholders with your actual data)
lambda_r = LamdaR(jj); % Your lambda_r value
a_prime = aprime(jj); % Your a_prime value

% Phi function interpolated over theta_values_deg in degrees
phi_vector_deg = Velocity_Direction_corrected(jj,1:360); % Your phi data as a vector, corresponding to each degree
theta_values_deg = 0:359; % A vector from 0 to 359 degrees
phi_function_deg = @(theta_deg) interp1(theta_values_deg, phi_vector_deg, theta_deg, 'spline');

% Alpha, Cl, and Cd data (replace these placeholders with your actual data)
alpha_data_deg = alpha_ext; % Your alpha values for Cd and Cl data, in degrees
Cd_data = CD_ext; % Your Cd data as a function of alpha, in degrees
Cl_data = CL_ext; % Your Cl data as a function of alpha, in degrees

% Define the range of theta_p values over which you want to plot
theta_p_deg = linspace(0,70,100); % for example from 0 to 180 degrees

% Call the function to plot the objective function
[theta_p_range_deg,objective_values] = plotObjectiveFunction(phi_vector_deg, lambda_r, a_prime, a(jj), alpha_data_deg, Cd_data, Cl_data, theta_p_deg);

[val,idx] = max(objective_values);
best_pitch(jj) = theta_p_range_deg(idx);


Elemental_CP(jj) = plotObjectiveFunctionSingleThetaP(phi_vector_deg, lambda_r, a_prime, a(jj), alpha_data_deg, Cd_data, Cl_data, best_pitch(jj));
QBlade_Elemental_CP(jj) = plotObjectiveFunctionSingleThetaP(phi_vector_deg, lambda_r, a_prime, a(jj), alpha_data_deg, Cd_data, Cl_data, Qblade_pitch(jj));

end

function [theta_p_range_deg,objective_values] = plotObjectiveFunction(phi_vector_deg, lambda_r, a_prime, a, alpha_data_deg, Cd_data, Cl_data, theta_p_range_deg)
   % Create phi interpolation function
    theta_deg_vec = 0:359; % Assuming phi_vector_deg is given for each degree
    phi_function_deg = @(theta_deg) interp1(theta_deg_vec, phi_vector_deg, theta_deg, 'spline');
  
   % Initialize the objective function value array
   objective_values = zeros(size(theta_p_range_deg));
  
   % Calculate the objective function value for each theta_p
   for i = 1:length(theta_p_range_deg)
       theta_p_deg = theta_p_range_deg(i);
       sum_term = 0;
       % Loop through all theta values where phi is defined
       for theta = theta_deg_vec
           % Calculate alpha for the current theta and theta_p
           alpha_deg = phi_function_deg(theta) - theta_p_deg;
          
           % Interpolate the values of Cd and Cl at the current alpha
           Cd_val = interp1(alpha_data_deg, Cd_data, alpha_deg, 'spline', 'extrap');
           Cl_val = interp1(alpha_data_deg, Cl_data, alpha_deg, 'spline', 'extrap');
          
           % Calculate the term of the objective function for the current theta
           if abs(phi_function_deg(theta))>2 && Cl_val > 0.05
               % Calculate the term of the objective function for the current theta
               temp = lambda_r^3 * a_prime * (1 - a) * (1 - ((Cd_val / Cl_val) * cotd(phi_function_deg(theta))));
               if abs(temp)>1000
                   term = 0;
               else
                   term = temp;

               end
           else
               term = 0;
           end
           sum_term = sum_term + term;
       end
      
       % Store the calculated value
       objective_values(i) = sum_term;
   end
  
   % Plot the objective function values against the theta_p range
   figure;
   plot(theta_p_range_deg, objective_values, 'b-', 'LineWidth', 2);
   xlabel('Pitch Angle');
   ylabel('Objective Function Value');
   title('Objective Function Variation with Pitch Angle');
   grid on;
end

function Elemental_CP = plotObjectiveFunctionSingleThetaP(phi_vector_deg, lambda_r, a_prime, a, alpha_data_deg, Cd_data, Cl_data, theta_p_deg)
    % Create phi interpolation function
    theta_deg_vec = 0:359; % Assuming phi_vector_deg is given for each degree
    phi_function_deg = @(theta_deg) interp1(theta_deg_vec, phi_vector_deg, theta_deg, 'spline');

    % Initialize the objective function value array for theta_deg_vec
    objective_values = zeros(size(theta_deg_vec));
    
    % Loop through all theta values where phi is defined
    for i = 1:length(theta_deg_vec)
        theta_deg = theta_deg_vec(i);
        % Calculate alpha for the current theta and the single theta_p
        alpha_deg = phi_function_deg(theta_deg) - theta_p_deg;

        % Interpolate the values of Cd and Cl at the current alpha
        Cd_val = interp1(alpha_data_deg, Cd_data, alpha_deg, 'spline', 'extrap');
        Cl_val = interp1(alpha_data_deg, Cl_data, alpha_deg, 'spline', 'extrap');

        % Calculate the term of the objective function for the current theta
        if abs(phi_function_deg(theta_deg)) > 2
            temp = lambda_r^3 * a_prime * (1 - a) * (1 - ((Cd_val / Cl_val) * cotd(phi_function_deg(theta_deg))));
            if abs(temp) > 1000
                term = 0;
            else
                term = temp;
            end
        else
            term = 0;
        end
        objective_values(i) = term; % Store the calculated value for the current theta
    end

    Elemental_CP = sum(objective_values)/360;
end



