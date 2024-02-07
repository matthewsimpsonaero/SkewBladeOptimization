function GenerateChordPlot(roverR,coverR_wake)
fig = figure();
fig.Position = [100 100 740 600];
plot(roverR,coverR_wake)
xlabel('r/R')
ylabel('c/R')
grid on
end

