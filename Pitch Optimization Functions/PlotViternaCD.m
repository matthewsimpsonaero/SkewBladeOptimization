function PlotViternaCD(alpha_ext,CD_ext,AOA_Curve, Drag_Curve)
fig = figure();
fig.Position = [100 100 740 600];
gg = plot(alpha_ext,CD_ext,'r',LineWidth=2,DisplayName='Viterna Extrapolation');
hold on
uu = plot(AOA_Curve, Drag_Curve,'m',LineWidth=2,DisplayName='Xfoil Output');
grid on 
grid(gca,'minor')
xlabel('Angle of Attack α (deg)',FontSize=16)
ylabel('CD',FontSize=16)
title('Viterna Expansion of Airfoil Curve \alpha \in [-180^\circ, 180^\circ]',FontSize=16)
xlim([-180,180])
ax = gca;
ax.XTick = -180:45:180;
ax.XTickLabel = {-180, -135, -90, -45, 0, 45, 90, 135, 180};
xline(0,'k-')
yline(0,'k-')
legend([uu,gg],fontsize=12)
ylim([0 1.4])
end

