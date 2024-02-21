function GenerateChordPlot(roverR,coverR_wake,coverR_no_wake)
fig = figure();
fig.Position = [100 100 740 600];
plot(roverR,coverR_wake,'b-',linewidth=2)
hold on
plot(roverR,coverR_no_wake,'r-',linewidth=2)
xlabel('r/R')
ylabel('c/R')
grid on
title('Nondimensionalized Chord vs. Nondimensionalized Span')
legend("With Wake Rotation",'Without Wake Rotation')
end

