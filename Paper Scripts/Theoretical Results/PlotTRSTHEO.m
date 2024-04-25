dataT20 = readmatrix("THAT20.xlsx");
dataT10 = readmatrix("THAT10.xlsx");
dataT0 = readmatrix("THAT0.xlsx");
dataQ20 = readmatrix("Q20.xlsx");
dataQ10 = readmatrix("Q10.xlsx");
dataQ0 = readmatrix("Q0.xlsx");


plot(dataT20(:,1),dataT20(:,2),'o-',LineWidth=1.5,Color=[0.4660 0.6740 0.1880],MarkerSize=9)
hold on
plot(dataT10(:,1),dataT10(:,2),"+-",LineWidth=1.5,Color=[0.4660 0.6740 0.1880],MarkerSize=9)
plot(dataT0(:,1),dataT0(:,2),"*-",LineWidth=1.5,Color=[0.4660 0.6740 0.1880],MarkerSize=9)
plot(dataQ20(:,1),dataQ20(:,2),"diamond-b",LineWidth=1.5,Color=[0 0.4470 0.7410],MarkerSize=9)
plot(dataQ10(:,1),dataQ10(:,2),"x-b",LineWidth=1.5,Color=[0 0.4470 0.7410],MarkerSize=9)
plot(dataQ0(:,1),dataQ0(:,2),"square-b",LineWidth=1.5,Color=[0 0.4470 0.7410],MarkerSize=9)
grid on
grid minor
title('Theoretical CP vs. TSR for Betz & QBlade Turbine',FontSize=14)
xlabel('Tip Speed Ratio (TSR)',FontSize=18)
ylabel('Power Coefficent (C_P)',FontSize=18)
legend('Novel \psi=20°','Novel \psi=10°','Novel \psi=0°','Betz \psi=20°','Betz \psi=10°','Betz \psi=0°')
ylim([0,.52])
