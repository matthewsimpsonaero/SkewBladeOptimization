af = readmatrix('airfoil.txt');

hold on
xlim([-.6 1])
ylim([-1 1])
grid on
grid minor
axis square

quiver(-0.5, 0, 0.54, 0, 'r', 'LineWidth', 2,'MaxHeadSize', .5) % X-axis vector
quiver(-0.5, -0.5, 0, 0.54, 'b', 'LineWidth', 2,'MaxHeadSize', .5) % Y-axis vector
quiver(-0.5, -0.5, .54, .54, 'g', 'LineWidth', 2,'MaxHeadSize', .5) % Y-axis vector

text(-0.1, 0.1, 'X', 'FontSize', 12, 'HorizontalAlignment', 'center')
text(-0.5, -0.55, 'Y', 'FontSize', 12,'HorizontalAlignment', 'center')

plot(af(:,1), af(:,2), 'k-', 'LineWidth', 2, 'HandleVisibility', 'off');
legend('X-Axis','Y-Axis','Resultant Vector')

th = linspace( 180, 225, 100);
R = .2;  %or whatever radius you want
x = R*cosd(th) ;
y = R*sind(th) ;
plot(x,y,'k-','HandleVisibility', 'off'); 

text(-0.25, -0.1, 'Φ', 'FontSize', 12,'HorizontalAlignment', 'center')
set(gca,'Yticklabel',[]) 
set(gca,'Xticklabel',[])
