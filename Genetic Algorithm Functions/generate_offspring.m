function Population = generate_offspring(Parent_Population,num_offspring)
% This function generates offspring by averaging the trails of the parents
% that are selected randomly. There is no mutation currently, but that
% could be a future feature.

Population = {}; % generate a new population variable
for i = 1:num_offspring % iterate until all of the offspring have been generated

    Parent1 = char(Parent_Population{randi([1,length(Parent_Population)])}); % select a parent 1 from the parent population
    Parent2 = char(Parent_Population{randi([1,length(Parent_Population)])}); % select a parent 2 from the parent population

    Parent1_camber = Parent1(1); % parse the traits of the parents
    Parent2_camber = Parent2(1);

    Parent1_camber_pos = Parent1(2);
    Parent2_camber_pos = Parent2(2);

    Parent1_thickness = Parent1(3:4);
    Parent2_thickness = Parent2(3:4);

    Offspring_camber = num2str(round(mean([str2double(Parent1_camber),str2double(Parent2_camber)]))); % average the camber of the parents
    Offspring_camber_pos = num2str(round(mean([str2double(Parent1_camber_pos),str2double(Parent2_camber_pos)]))); % average the camber position of the parents
    Offspring_thickness = round(mean([str2double(Parent1_thickness),str2double(Parent2_thickness)])); % average the thickness of the parents
    if Offspring_thickness < 10
        % if the thickness is less than 10, then a leading 0 needs to be
        % added. For example, for a camber position of 5, the output is
        % required to be 05 in order to give each NACA airfoil 4 digits
        Offspring_thickness = [num2str(0),num2str(Offspring_thickness)];
    else
        Offspring_thickness = num2str(Offspring_thickness);
    end
    % assemble the offspring airfoil code. This will look like 1230, 0010, 2604 etc.
    Offspring = [Offspring_camber,Offspring_camber_pos,Offspring_thickness];
    
    % store the new airfoil code with the others 
    Population{i} = Offspring;
end

    

end

