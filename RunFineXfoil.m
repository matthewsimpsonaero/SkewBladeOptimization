function [Maxlift,bestAOA] = RunFineXfoil(airfoil,Element_Reynolds,Element_Mach)

    fine_xfoil = 50;
    points_fine = linspace(0,20,fine_xfoil);
    [Polar] = Airfoil_Runner(airfoil,Element_Reynolds,Element_Mach,points_fine);
    [Maxlift,idx] = max(Polar.CL);
    bestAOA = Polar.Alpha(idx);

end

