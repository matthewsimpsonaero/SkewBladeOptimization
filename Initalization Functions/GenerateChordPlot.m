function GenerateChordPlot(roverR,coverR_wake,coverR_no_wake)
fig = figure();
fig.Position = [100 100 740 600];
plot(roverR,coverR_wake,'b-',linewidth=2)
hold on
plot(roverR,coverR_no_wake,'r-',linewidth=2)
xlabel('Non dimensional Radius (r/R)',fontsize = 16)
ylabel('Non dimensional Chord Length (c/R)',fontsize = 16)
grid on
grid(gca,'minor')

title('Nondimensionalized Chord vs. Nondimensionalized Span',fontsize = 16)
legend("With Wake Rotation",'Without Wake Rotation',fontsize = 16)
end

