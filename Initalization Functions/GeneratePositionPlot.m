function GeneratePositionPlot(Position,n)
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
end

