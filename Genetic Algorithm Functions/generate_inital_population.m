function Airfoils = generate_inital_population(Population_Size)
% This function generates the inital population by randomly selecting a
% camber (0-9), a camber position (0-9), and a thickness (0-40). The
% airfoils are stored in a cell array called 'Airfoils'

Airfoils = {}; % initalize cell array

for i = 1:Population_Size % iterate until the population size has been reached
    Camber = randi([0,9]); % select camber
    if Camber>0
        Camber_Position = randi([0 9]); % if there is camber, select a camber position
    else
        Camber_Position = 0; % if there is no camber, then the camber position is irrelevant
    end

    Thickness = randi([10 40]); % select a thickness
    if Thickness < 10
        % if the thickness is less than 10, then a leading 0 needs to be
        % added. For example, for a camber position of 5, the output is
        % required to be 05 in order to give each NACA airfoil 4 digits
        Thickness = [num2str(0),num2str(Thickness)]; 
    end
    % assemble the airfoil code. This will look like 1230, 0010, 2604 etc.
    Af = [num2str(Camber),num2str(Camber_Position),num2str(Thickness)];
    
    % store the airfoil code with the others 
    Airfoils{i} = Af;
end

