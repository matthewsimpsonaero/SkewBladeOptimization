% % compute solidarity
% R_values = Element_location+hub_radius;
% R_total = Blade_length+hub_radius;
% integral_C = trapz(R_values, coverR_wake*(R_total));
% sigma = ((integral_C*blade_number) +(pi*hub_radius^2)) / (pi * R_total^2); %Wiley 3.107

% Compute axial and tangential induction factors

% for i = 1:n+1
% aprime(i) = 1/(((4*cosd(Average_Velocity_Direction(i)))/(sigma*GuessCL))-1); % Wiley 3.89
% a(i) = -((aprime(i)/((sigma*GuessCL*cosd(Average_Velocity_Direction(i)))/(4*LamdaR(i)*sind(Average_Velocity_Direction(i)))))-1); %wiley 3-83
% end
% 
% % Compute Velocity using the induction factors
% 
% for i = 1:n+1
% V_parallel_corrected(i,:) = (V_parallel(i,:).*(1+aprime(i)));
% V_perpindicular_corrected(i,:) = (V_perpindicular(i,:).*(1-a(i)));
% Velocity_Direction_corrected(i,1:360) = atan2d(V_perpindicular_corrected(i,:), V_parallel_corrected(i,:));
% end
% 
% if generate_plots
% fig = figure();
% fig.Position = [100 100 740 600];
% for i = 1:n+1
% plot(1:360,sqrt(V_parallel_corrected(i,:).^2+V_perpindicular_corrected(i,:) .^2))
% hold on
% end
% grid on 
% grid(gca,'minor')
% xlabel('° (deg)')
% ylabel('Seen Velocity (m/s)')
% title('Seen Velocity vs. Rotation Angle for Skewed Turbine')
% 
% fig = figure();
% fig.Position = [100 100 740 600];
% for i = 1:n+1
% plot(1:360,Velocity_Direction_corrected(i,:))
% hold on
% end
% grid on 
% grid(gca,'minor')
% xlabel('° (deg)')
% ylabel('Velocity Direction (°)')
% title('Velocity Direction vs. Rotation Angle for Skewed Turbine')
% end