function [theta_p_range_deg,objective_values] = plotObjectiveFunction(phi_vector_deg, lambda_r, a_prime, a, alpha_data_deg, Cd_data, Cl_data, theta_p_range_deg)
    % Create phi interpolation function
    theta_deg_vec = 0:359; % Assuming phi_vector_deg is given for each degree
    phi_function_deg = @(theta_deg) interp1(theta_deg_vec, phi_vector_deg, theta_deg, 'spline');
  
   % Initialize the objective function value array
   objective_values = zeros(size(theta_p_range_deg));
  
   % Calculate the objective function value for each pitch angle
   for i = 1:length(theta_p_range_deg)
       theta_p_deg = theta_p_range_deg(i); % pitch the pitch value
       sum_across_azimtuh = 0; % initalize the summation across the azimuth sweep
 
       for theta = theta_deg_vec % iterate through each azimuth angle
           % Calculate the angle of attack at this azimuth angle
           alpha_deg = phi_function_deg(theta) - theta_p_deg;
          
           % Calcualte the CL and CD at this azimuth angle
           Cd_val = interp1(alpha_data_deg, Cd_data, alpha_deg, 'spline', 'extrap');
           Cl_val = interp1(alpha_data_deg, Cl_data, alpha_deg, 'spline', 'extrap');
          
           % if phi is close to zero (asymptote) or CL is close to zero
           % (asymptote), skip this azimuth angle in the calculation
           if abs(phi_function_deg(theta))>2 && Cl_val > 0.05
               % Use our pitch optimization formula from the AIAA paper
               temp = lambda_r^3 * a_prime(theta+1) * (1 - a(theta+1)) * (1 - ((Cd_val / Cl_val) * cotd(phi_function_deg(theta))));
               if abs(temp)>1000 % if the result is super large (bad) skip this azimuth angle in the calculation
                   term = 0;
               else
                   term = temp; % if no issues, store the value
               end
           else
               term = 0;
           end
           sum_across_azimtuh = sum_across_azimtuh + term; % add this contribution to the CP to the previous summation
       end
      
       % Store the total value
       objective_values(i) = sum_across_azimtuh;
   end
  
   % Plot the objective function values against the candidate pitch range
   fig = figure();
   fig.Position = [100 100 740 600];
   plot(theta_p_range_deg, objective_values, 'b-', 'LineWidth', 2);
   hold on
   [value,index] = max(objective_values);
   plot(theta_p_range_deg(index), value, 'r.', 'MarkerSize',30);
   xlabel('Element Pitch Angle (deg)',FontSize=16);
   ylabel('Objective Function Value',FontSize=16);
   title('Objective Function Variation with Pitch Angle',FontSize=16);
   grid on;
   legend('Objective Function','Best Pitch Location',FontSize=16)
end
