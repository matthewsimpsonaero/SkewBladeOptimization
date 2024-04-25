% Matthew Simpson
% Pitch Optimizer using QBlade Outputs
% 2-21-2024
clc;clear, close all

%% Code Setup

% This process is extremely manual. An outline of how to run this code will
% be provided. This process could be automated, but time time complexity of
% that may outweight the time saving of doing so in the short term.

addpath(genpath('./Pitch Optimization Functions')) 

% Add the path to the QBlade Output files
addpath(genpath('./Twist Optimization 20 Deg')) 

Element_Num = 5; % modify this to select which element that is being run
Airfoil_Code = '5510'; % modify this to select the elements airfoil
Average_Mach = 0.0797; % modify this to match the average mach number of the element
lambda_r = 5.23; % modify this to match the lamda_r (elemental TSR) of the element
QBlade_Pitch_Value = -0.194; % modify this to match the pitch value supplied from Qblade for this run
theta_p_deg = linspace(-20,70,200); % Define the range of theta_p values over which you want to run the optimiation over

%% Pull the reynolds number from the file

Reynolds_Distrobution = readmatrix(sprintf('/Element %d/Reynolds_%d.txt',Element_Num,Element_Num)); % get the reynolds vs. azimuth angle from the file
Average_Reynolds_Number = mean(Reynolds_Distrobution(361:720,2)); % find the average reynolds number from the file


%% Run Xfoil to get the CL and CD vs. Angle of attack for the airfoil

% Run fine xfoil to get a high fidelity CL and CD curve
[Maxlift,bestAOA,Lift_curve,Drag_curve,AOA_curve] = RunFineXfoil(Airfoil_Code,Average_Reynolds_Number(1),Average_Mach); %$$$$$$$$$$$$ X2

% store the output to their respective variables
Lift_Curve = Lift_curve;
Drag_Curve = Drag_curve;
AOA_Curve = AOA_curve;


%% Perform Vitera Extrapolation 

[alpha_ext, CL_ext, CD_ext] = viterna_extrapolation(AOA_Curve, Lift_curve, Drag_Curve); %perform vinerna expansion 

% plot the extrapolations
PlotViternaCD(alpha_ext,CD_ext,AOA_Curve, Drag_Curve)
PlotViternaCL(alpha_ext,CL_ext,AOA_Curve, Lift_Curve)


%% Setup the Pitch Optimization

% Pull the tangential induction factor curve from the QBlade Data
Tangential_Induction_Curve = readmatrix(sprintf('/Element %d/Tangential_Induction_%d.txt',Element_Num,Element_Num));
a_prime = Tangential_Induction_Curve(361:720,2); 

% Pull the axial induction factor curve from the QBlade Data
Axial_Induction_Curve = readmatrix(sprintf('/Element %d/Axial_Induction_%d.txt',Element_Num,Element_Num)); 
a  = Axial_Induction_Curve(361:720,2);

% Pull the angle of attack curve from the QBlade Data
angleattack = readmatrix( sprintf('/Element %d/AOA_%d.txt',Element_Num,Element_Num));

% The phi vector is the angle of attack + the QBlade pitch value determined
% through the Qblade optimization
phi_vector_deg = angleattack(361:720,2)+QBlade_Pitch_Value;

% Call the optimzation function
[theta_p_range_deg,objective_values] = plotObjectiveFunction(phi_vector_deg, lambda_r, a_prime, a, alpha_ext, CD_ext, CL_ext, theta_p_deg);

% find the maximum value from the objective function, and determine the
% pitch angle that generated the maximum values
[val,idx] = max(objective_values);
best_pitch = theta_p_range_deg(idx);

fprintf('The best pitch is: %0.3f°\n',best_pitch)
