clc;clear

addpath(genpath('./QBlade Pitch Files Final')) 
Reynolds_1 = readmatrix('/Element 10/Reynolds_10.txt'); %%$$$$$$$$$$$$$
Average_Reynolds_Number = mean(Reynolds_1(1:359,2));

plot(Reynolds_1(1:359,1),Reynolds_1(1:359,2))

%% Load best airfoils

fileID = fopen('Emelia Results/integration/results_final.txt', 'r');
fgetl(fileID);
C = textscan(fileID, '%s', 'Delimiter', '\n');
dataCellArray = C{1};
fclose(fileID);

for i = 1:10
    temp = dataCellArray{i};
    if i<10
        NACAnum(i) = str2double(temp(14:17));
    else
        NACAnum(i) = str2double(temp(15:18));
    end
end

liftofElement = {};
AOAofElement ={};
Lift_curves = {};
Drag_curves = {};
AOA_curves = {};

for i = 1:1
[Maxlift,bestAOA,Lift_curve,Drag_curve,AOA_curve] = RunFineXfoil('6710',Average_Reynolds_Number(1),0.082); %$$$$$$$$$$$$ X2

liftofElement{i} = Maxlift;
AOAofElement{i} = bestAOA;
Lift_curves{i} = Lift_curve;
Drag_curves{i} = Drag_curve;
AOA_curves{i} = AOA_curve;
end


%% Determine twist

for jj = 1:1
[alpha_ext, CL_ext, CD_ext] = viterna_extrapolation(AOA_curves{jj}, Lift_curves{jj}, Drag_curves{jj}); %perform vinerna expansion 

% Define the known values (replace these placeholders with your actual data)
lambda_r = 5.23; % Your lambda_r value $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
t1 = readmatrix('/Element 10/Tangential_Induction_10.txt'); %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
a_prime = t1(1:360,2); % Your a_prime value
t2 = readmatrix('/Element 10/Axial_Induction_10.txt'); %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
a  = t2(1:360,2);

% Phi function interpolated over theta_values_deg in degrees
angleattack = readmatrix('/Element 10/AOA_10.txt'); %$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
phi_vector_deg = angleattack(1:360,2)+3.264; % $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
theta_values_deg = 0:359; % A vector from 0 to 359 degrees
phi_function_deg = @(theta_deg) interp1(theta_values_deg, phi_vector_deg, theta_deg, 'spline');

% Alpha, Cl, and Cd data (replace these placeholders with your actual data)
alpha_data_deg = alpha_ext; % Your alpha values for Cd and Cl data, in degrees
Cd_data = CD_ext; % Your Cd data as a function of alpha, in degrees
Cl_data = CL_ext; % Your Cl data as a function of alpha, in degrees

% Define the range of theta_p values over which you want to plot
theta_p_deg = linspace(0,70,200); % for example from 0 to 180 degrees

% Call the function to plot the objective function
[theta_p_range_deg,objective_values] = plotObjectiveFunction(phi_vector_deg, lambda_r, a_prime, a, alpha_data_deg, Cd_data, Cl_data, theta_p_deg);

[val,idx] = max(objective_values);
best_pitch(jj) = theta_p_range_deg(idx);


%Elemental_CP(jj) = plotObjectiveFunctionSingleThetaP(phi_vector_deg, lambda_r, a_prime, a, alpha_data_deg, Cd_data, Cl_data, best_pitch(jj));
%QBlade_Elemental_CP(jj) = plotObjectiveFunctionSingleThetaP(phi_vector_deg, lambda_r, a_prime, a, alpha_data_deg, Cd_data, Cl_data, Qblade_pitch(jj));

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
               temp = lambda_r^3 * a_prime(theta+1) * (1 - a(theta+1)) * (1 - ((Cd_val / Cl_val) * cotd(phi_function_deg(theta))));
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

% function Elemental_CP = plotObjectiveFunctionSingleThetaP(phi_vector_deg, lambda_r, a_prime, a, alpha_data_deg, Cd_data, Cl_data, theta_p_deg)
%     % Create phi interpolation function
%     theta_deg_vec = 0:359; % Assuming phi_vector_deg is given for each degree
%     phi_function_deg = @(theta_deg) interp1(theta_deg_vec, phi_vector_deg, theta_deg, 'spline');
% 
%     % Initialize the objective function value array for theta_deg_vec
%     objective_values = zeros(size(theta_deg_vec));
%     
%     % Loop through all theta values where phi is defined
%     for i = 1:length(theta_deg_vec)
%         theta_deg = theta_deg_vec(i);
%         % Calculate alpha for the current theta and the single theta_p
%         alpha_deg = phi_function_deg(theta_deg) - theta_p_deg;
% 
%         % Interpolate the values of Cd and Cl at the current alpha
%         Cd_val = interp1(alpha_data_deg, Cd_data, alpha_deg, 'spline', 'extrap');
%         Cl_val = interp1(alpha_data_deg, Cl_data, alpha_deg, 'spline', 'extrap');
% 
%         % Calculate the term of the objective function for the current theta
%         if abs(phi_function_deg(theta_deg)) > 2
%             temp = lambda_r^3 * a_prime(theta_deg+1) * (1 - a(theta_deg+1)) * (1 - ((Cd_val / Cl_val) * cotd(phi_function_deg(theta_deg))));
%             if abs(temp) > 1000
%                 term = 0;
%             else
%                 term = temp;
%             end
%         else
%             term = 0;
%         end
%         objective_values(i) = term; % Store the calculated value for the current theta
%     end
% 
%     Elemental_CP = sum(objective_values)/360;
% end