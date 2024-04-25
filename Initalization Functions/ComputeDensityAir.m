function Density_air = ComputeDensityAir(air_temp,atmos_pressure,R,molar_mass_air)
Density_air = atmos_pressure/((R*1000/molar_mass_air)*(air_temp+273.15)); % Ideal gas law
end

