function GenerateSkewVelocityPlot(Skew_Velocity,n)
fig = figure();
fig.Position = [100 100 740 600];
for i = 1:n+1
    plot(1:360,Skew_Velocity(i,:),LineWidth=2,DisplayName=sprintf('E. %d',i-1))
    hold on
end
grid on
grid(gca,'minor')
xlabel('Azimuth Angle θ (deg)',FontSize=16)
ylabel('Skew Velocity (m/s)',FontSize=16)
title('Skew Velocity vs. Rotation Angle for Skewed Turbine',FontSize=16)
leg = legend(Location="southoutside",Orientation='horizontal');
leg.ItemTokenSize = [10,6];
xlim([0,360])
ax = gca;
ax.XTick = 0:45:360;
ax.XTickLabel = {0, 45, 90, 135, 180, 225, 270, 315, 360};
end

