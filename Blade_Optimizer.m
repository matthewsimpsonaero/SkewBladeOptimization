% Matthew Simpson
% Turbine Blade Optimizer in Skew
% 1-27-2024
clc;clear

%% Input parameters

V_inf = 5; %m/s
Blade_length = 0.15; %m
skew_angle = 30; %deg
blade_number = 3;
TSR = (4*pi/blade_number)*1.25;  %https://users.wpi.edu/~cfurlong/me3320/DProject/Ragheb_OptTipSpeedRatio2014.pdf
Chord_inital_guess = 0.06; %cm
n = 10; %number of blade element segments
hub_radius = 0.0127; %meters

%% Flow Paramters

air_temp = 20; %celcius
atmos_pressure = 101325; %Pa
R = 8.314459; %j/mol*K
molar_mass_air = 28.9628; %g/mol
dynamic_viscosity_air = 1.8205e-5; % kg/m*s
gamma_air = 1.4;

Density_air = atmos_pressure/((R*1000/molar_mass_air)*(air_temp+273.15));

%% Determination of reynolds number
close all

Blade_element_index = 1:n;
Element_location = linspace(hub_radius,Blade_length+hub_radius,n);
omega = (TSR*V_inf)/ (Blade_length+hub_radius);

for i = 1:n
for theta = 1:360
    Position(i,theta) = Element_location(i)*sind(skew_angle)*cosd(theta);
    Skew_Velocity(i,theta) = -omega*Element_location(i)*sind(skew_angle)*sind(theta);
end 
end

fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,Position(i,:))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Axial Position (m)')
title('Axial Position vs. Rotation Angle for Skewed Turbine')

fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,Skew_Velocity(i,:))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Skew Velocity (m/s)')
title('Skew Velocity vs. Rotation Angle for Skewed Turbine')

V_total_static = Skew_Velocity+(V_inf*cosd(skew_angle));

fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,V_total_static(i,:))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Total Axial Velocity (m/s)')
title('Total Axial Velocity vs. Rotation Angle for Skewed Turbine')

theta = 1:360;
for i = 1:n
V_parallel(i,1:360) = V_inf*cosd(skew_angle)-omega*Element_location(i)*sind(skew_angle)*sind(theta);
V_perpindicular(i,1:360) = V_inf*sind(skew_angle)*sind(theta) + omega*Element_location(i);
Velocity_Direction(i,1:360) = atan2d( V_perpindicular(i,1:360), V_parallel(i,1:360));
end

fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n
plot(1:360,sqrt(V_parallel(i,:).^2+V_perpindicular(i,:).^2))
hold on
end
grid on 
grid(gca,'minor')
xlabel('° (deg)')
ylabel('Seen Velocity (m/s)')
title('Seen Velocity vs. Rotation Angle for Skewed Turbine')

for i = 1:n
Reynolds_Number(i,1:360) = (Density_air*sqrt(V_parallel(i,:).^2+V_perpindicular(i,:).^2)*Chord_inital_guess)/dynamic_viscosity_air;
end

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

Average_Reynolds_Number = mean(Reynolds_Number');

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

Average_Velocity_Direction = mean(Velocity_Direction');

for i = 1:n
Direction_range(i) = max(Velocity_Direction(i,:))-min(Velocity_Direction(i,:));
end

for i = 1:n
Average_Velocity(i) = mean(sqrt(V_parallel(i,:).^2+V_perpindicular(i,:).^2));
end

for i = 1:n
Average_Mach(i) = mean(sqrt(V_parallel(i,:).^2+V_perpindicular(i,:).^2) ./ sqrt(gamma_air*(R*1000/molar_mass_air)*(air_temp+273.15)));
end

%% Xfoil Selection Using a Genetic Algorithm

fid2 = fopen( 'results.txt', 'wt' );
fid3 = fopen( 'Generational_Best.txt', 'wt' );
fprintf(fid2, 'Element Number\t\tAirfoil\t\tCL\t\tAOA\n');
fprintf(fid3, 'Element Number\t\tAirfoil\t\tAOA\n');

% course Xfoil
%[~,txt] = xlsread('NACA Candidates.xlsx');

% starting_population = 3240/5; % 20 percent of the total subset
% num_generations = 10;
% num_parents = starting_population/2;
% num_offspring = starting_population/2;

starting_population = 40; % 20 percent of the total subset
num_generations = 2;
num_parents = starting_population/2;
num_offspring = starting_population/2;


for element_num = 1:n

    Element_Reynolds = Average_Reynolds_Number(element_num);
    Element_Mach = Average_Mach(element_num);

    Population = generate_inital_population(starting_population);
    for gen = 1:num_generations
        for i = 1:length(Population)
            NACA_Number = Population{i};
            course_xfoil = 25; % number of points
            points_course = linspace(0,20,course_xfoil);
            try
                fprintf('G%d:(%d/%d)Running NACA%s\n',gen,i,length(Population),NACA_Number)
                [Polar] = Airfoil_Runner(NACA_Number,Element_Reynolds,Element_Mach,points_course);
                Lift(i) = max(Polar.CL);
            catch
                Lift(i) = 0;
            end
        end
        [sorted_values, sorted_indexes] = sort(Lift, 'descend');
        top_values = sorted_values(1:num_parents);
        indexes = sorted_indexes(1:num_parents);
        Parents = {};
        for i = 1:num_parents
            Parents{i} = Population(indexes(i));
        end
        fprintf(fid3, '%d\t\tNACA%s\t\t%.3f\t\t%d\n',element_num,Population{sorted_indexes(1)},sorted_values(1),Element_Reynolds);  
        if gen~=num_generations
            Population = generate_offspring(Parents,num_offspring);
            clear Lift
        end
    end

    fine_xfoil = 50;
    points_fine = linspace(0,20,fine_xfoil);
    [Polar] = Airfoil_Runner(Population{indexes(1)},Element_Reynolds,Element_Mach,points_fine);
    [Maxlift,idx] = max(Polar.CL);
    bestAOA = Polar.Alpha(idx);
    
    
    fprintf('Best Airfoil: %s\n',Population{indexes(1)})
    fprintf(fid2, '%d\t\tNACA%s\t\t%.3f\t\t%.3f\t\t%d\n',element_num,Population{indexes(1)},Maxlift,bestAOA,Element_Reynolds);
end



