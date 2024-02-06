function Airfoils = generate_inital_population(Population_Size)

Airfoils = {};

for i = 1:Population_Size
    Camber = randi([0,9]);
    if Camber>0
        Camber_Position = randi([0 9]);
    else
        Camber_Position = 0;
    end
    Thickness = randi([1 40]);
    if Thickness < 10
        Thickness = [num2str(0),num2str(Thickness)];
    end
    Af = [num2str(Camber),num2str(Camber_Position),num2str(Thickness)];

    Airfoils{i} = Af;



end

