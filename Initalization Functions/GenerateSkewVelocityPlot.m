function GenerateSkewVelocityPlot(Skew_Velocity,n)
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
end

