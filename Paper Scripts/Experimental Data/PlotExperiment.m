close all

dataT20 = readmatrix("THAT20Output.csv");
dataT10 = readmatrix("THAT10Output.csv");
dataT0 = readmatrix("THAT0Output.csv");
dataQ20 = readmatrix("QB20Output.csv");
dataQ10 = readmatrix("QB10Output.csv");
dataQ0 = readmatrix("QB0Output.csv");


fig = figure();
fig.Position = [100 100 1540 600];
subplot(1,3,1)
plot(dataT20(:,7),dataT20(:,8),'-o',LineWidth=1.5,MarkerSize=10,Color=[0, 0.4470, 0.7410])
hold on
plot(dataQ20(:,7),dataQ20(:,8),'-diamond',LineWidth=1.5,MarkerSize=10,Color=[0.4940, 0.1840, 0.5560])
grid on
grid minor
ylim([0 0.04])
xlim([4.4,5.8])
legend('Novel \psi=20°','Betz \psi=20°',Orientation='horizontal',FontSize=14)

subplot(1,3,2)
plot(dataT10(:,7),dataT10(:,8),'-+',LineWidth=1.5,MarkerSize=10,Color=[0.6350, 0.0780, 0.1840])
hold on 
plot(dataQ10(:,7),dataQ10(:,8),'-x',LineWidth=1.5,MarkerSize=10,Color=[0.4660, 0.6740, 0.1880])
grid on
grid minor
ylim([0 0.04])
xlim([4.4,5.8])
legend('Novel \psi=10°','Betz \psi=10°',Orientation='horizontal',FontSize=14)

subplot(1,3,3)
plot(dataT0(:,7),dataT0(:,8),'-*',LineWidth=1.5,Color=[0.9290, 0.6940, 0.1250],MarkerSize=10)
hold on
plot(dataQ0(:,7),dataQ0(:,8),'-s',LineWidth=1.5,MarkerSize=10,Color=[0.3010, 0.7450, 0.9330])
grid on
grid minor
ylim([0 0.04])
xlim([4.4,5.8])
legend('Novel \psi=0°','Betz \psi=0°',Orientation='horizontal',FontSize=14)

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
lh = ylabel(han,'Power Coefficent (C_P)',FontSize=18);
lh.Position(1) = lh.Position(1)-0.02;
xlabel(han,'Tip Speed Ratio (TSR)',FontSize=18)
lt = title(han,'Experimental CP vs. TSR for Betz & Novel Turbine',FontSize=18)
lt.Position(2) = lt.Position(2)+0.02;













%yline(max(dataQ10(:,8)),'m--',LineWidth=1.5)
%yline(max(dataT10(:,8)),'--',LineWidth=1.5,Color=[0.8500 0.3250 0.0980])
grid on
grid minor
% title('Experimental CP vs. TSR for Betz & QBlade Turbine',FontSize=14)
xlabel('Tip Speed Ratio (TSR)',FontSize=18)
% ylabel('Power Coefficent (C_P)',FontSize=18)
% legend('Novel \psi=20°','Novel \psi=10°','Novel \psi=0°','Betz \psi=20°','Betz \psi=10°','Betz \psi=0°','Maximum Betz Power','Maximum Novel Power',Location='southeast')
% xlim([4.4,6])
