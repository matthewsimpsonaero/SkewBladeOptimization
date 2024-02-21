function [Maxlift,bestAOA,Lift_curve,Drag_curve,AOA_curve,maxLoverDint] = RunFineXfoil(airfoil,Element_Reynolds,Element_Mach)

    fine_xfoil = 70;
    points_fine = linspace(-15,20,fine_xfoil);
    fprintf('Running NACA%s\n',airfoil)
    [Polar] = Airfoil_Runner(airfoil,Element_Reynolds,Element_Mach,points_fine);
    [Maxlift,idx] = max(Polar.CL);
    bestAOA = Polar.Alpha(idx);
    Lift_curve = Polar.CL;
    Drag_curve = Polar.CD;
    AOA_curve = Polar.Alpha;
    maxLoverDint =trapz(Polar.Alpha,(Polar.CL./Polar.CD));

end

