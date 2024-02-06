course_xfoil = 10; % number of points
    points_course = linspace(0,20,course_xfoil);
    
    for j = 1:length(txt)
        NACA_Number = txt{j};
        [Polar] = Airfoil_Runner(NACA_Number,test_re,test_mach,points_course);
        try
        lift(:,j) = Polar.CL;
        catch
        lift(:,j) = zeros(course_xfoil,1);
        end
    end
 
    % selecting 3 candidates for a higher order simulation 
    maximum_lift_course = max(lift);
    [sortedValues, sortedIndices] = sort(lift, 'descend');
    topThreeValues = sortedValues(1:3);
    topThreeLocations = sortedIndices(1:3);
    
    % Fine Xfoil
    
    fine_xfoil = 20; % number of points
    points_fine = linspace(0,20,fine_xfoil);
    
    for j = 1:length(topThreeLocations)
        Fine_NACA_Number = txt{topThreeLocations(j)};
        [Polar] = Airfoil_Runner(Fine_NACA_Number,test_re,test_mach,points_fine);
        try
        fine_lift(:,j) = Polar.CL;
        catch
        fine_lift(:,j) = zeros(length(fine_lift),1);
        end
    end
      
    
    maximum_lift_fine = max(fine_lift);
    [sortedValues, sortedIndices] = sort(maximum_lift_fine, 'descend');
    Top_val = sortedValues(1);
    Top_location = sortedIndices(1);
    best_airfoil = txt{Top_location};
    best_lift_sweep = fine_lift(:,Top_location);
    [maxValues, rowIndices] = max(best_lift_sweep);
    optimal_AOA = points_fine(rowIndices);