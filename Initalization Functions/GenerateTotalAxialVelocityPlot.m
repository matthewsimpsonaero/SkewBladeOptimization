function GenerateTotalAxialVelocityPlot(V_total_static,n)
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


end

