function Airfoils = generate_inital_population(Population_Size)

Airfoils = {};
counter = 1;
for i = 0:9
for j = 5:9
    Camber = i;
    if i == 0
    Camber_Position = 0;
    else
        Camber_Position = j;
    end
    Thickness = 10;
    Af = [num2str(Camber),num2str(Camber_Position),num2str(Thickness)];

    Airfoils{counter} = Af;
    counter = counter + 1;
end
end

% for i = 1:Population_Size
%     Camber = randi([0,9]);
%     if Camber>0
%         Camber_Position = randi([0 9]);
%     else
%         Camber_Position = 0;
%     end
%     Thickness = randi([10 40]);
%     if Thickness < 10
%         Thickness = [num2str(0),num2str(Thickness)];
%     end
%     Af = [num2str(Camber),num2str(Camber_Position),num2str(Thickness)];
% 
%     Airfoils{i} = Af;



end

