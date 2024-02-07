function Population = generate_offspring(Parent_Population,num_offspring)

Population = {};
for i = 1:num_offspring
    Parent1 = char(Parent_Population{randi([1,length(Parent_Population)])});
    Parent2 = char(Parent_Population{randi([1,length(Parent_Population)])});

    Parent1_camber = Parent1(1);
    Parent2_camber = Parent2(1);

    Parent1_camber_pos = Parent1(2);
    Parent2_camber_pos = Parent2(2);

    Parent1_thickness = Parent1(3:4);
    Parent2_thickness = Parent2(3:4);

    Offspring_camber = num2str(round(mean([str2double(Parent1_camber),str2double(Parent2_camber)])));
    Offspring_camber_pos = num2str(round(mean([str2double(Parent1_camber_pos),str2double(Parent2_camber_pos)])));
    Offspring_thickness = round(mean([str2double(Parent1_thickness),str2double(Parent2_thickness)]));
    if Offspring_thickness < 10
        Offspring_thickness = [num2str(0),num2str(Offspring_thickness)];
    else
        Offspring_thickness = num2str(Offspring_thickness);
    end

    Offspring = [Offspring_camber,Offspring_camber_pos,Offspring_thickness];
    Population{i} = Offspring;
end

    

end

