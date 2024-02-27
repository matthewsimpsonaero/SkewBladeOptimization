function [alpha_sweep,CL_export, CD_export] = viterna_extrapolation(alpha, CL, CD)
%https://www.researchgate.net/publication/317305962_Airfoil_Lift_and_Drag_Extrapolation_with_Viterna_and_Montgomerie_Methods
alpha_ext = 1:90;
AR = 3;
[Cl_stall,idx] = max(CL);
alpha_stall = alpha(idx);
Cd_stall = CD(idx);
CdMax = 1.11+0.018*AR;
A1 = CdMax/2;
B1 = CdMax;
A2 = (Cl_stall-CdMax.*sind(alpha_stall)*cosd(alpha_stall))*(sind(alpha_stall)/cosd(alpha_stall)^2);
B2 = (Cd_stall-CdMax*sind(alpha_stall)^2)/(cosd(alpha_stall));

for i = 1:length(alpha_ext)
    if alpha_ext(i) <= max(alpha)
        CL_ext(i) = interp1(alpha, CL, alpha_ext(i), 'linear', 'extrap');
        CD_ext(i) = interp1(alpha, CD, alpha_ext(i), 'linear', 'extrap');
    else 
        CL_ext(i) = A1.*sind(2*alpha_ext(i))+(A2*(cosd(alpha_ext(i)).^2./sind(alpha_ext(i))));
        CD_ext(i) = B1.*sind(alpha_ext(i)).^2 + B2.*cosd(alpha_ext(i));
    end
end
RHSCL = [CL_ext,-fliplr(CL_ext)];
RHSCD = [CD_ext,fliplr(CD_ext)];

TotalCL = [-fliplr(RHSCL),0,RHSCL];
TotalCD = [fliplr(RHSCD),0,RHSCD];
alpha_sweep = -180:180;

for i = 1:length(alpha_sweep)
    if alpha_sweep(i) >= min(alpha) && alpha_sweep(i) < max(alpha)
        CL_export(i) = interp1(alpha, CL, alpha_sweep(i), 'linear', 'extrap');
        CD_export(i) = interp1(alpha, CD, alpha_sweep(i), 'linear', 'extrap');
    else
         CL_export(i) = TotalCL(i);
         CD_export(i) = TotalCD(i);
    end

end
end
