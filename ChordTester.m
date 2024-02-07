TSR = 7;
GuessCL = 1;
blade_number=3;
n=10;
roverR =linspace(0+(1/n),1,10);
LamdaR = TSR*roverR;
wind_angle =[43.6,25.5,17.6,13.4,10.8,9,7.7,6.8,6,5.4];

for i = 1:n
chord(i) = (8*pi*roverR(i)*sind(wind_angle(i)))/(3*blade_number*GuessCL*LamdaR(i));
end
chord

plot(roverR,chord)