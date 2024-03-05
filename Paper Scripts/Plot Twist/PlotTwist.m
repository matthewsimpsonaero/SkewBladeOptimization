close all
data20 = readmatrix('Log Series 20 deg.xlsx');
data30 = readmatrix('Power Series 30 deg.xlsx');

figure()
plot(data20(:,1),data20(:,3),'bo-',LineWidth=2,MarkerSize=10)
hold on
plot(data20(:,1),data20(:,4),'g*-',LineWidth=2,MarkerSize=10)
grid on
grid minor
xlabel('Nondimensional Radius (r/R)',FontSize=16)
ylabel('Elemental Pitch (°)',FontSize=16)
legend('Novel Twist','Betz Twist')
title('20° Skew Optimized Blade Twist Function',FontSize=16)
ylim([-1,55])

figure()
plot(data30(:,1),data30(:,3),'bo-',LineWidth=2,MarkerSize=10)
hold on
plot(data30(:,1),data30(:,4),'g*-',LineWidth=2,MarkerSize=10)
grid on
grid minor
xlabel('Nondimensional Radius (r/R)',FontSize=16)
ylabel('Elemental Pitch (°)',FontSize=16)
title('30° Skew Optimized Blade Twist Function',FontSize=16)
legend('Novel Twist','Betz Twist')
ylim([-1,55])