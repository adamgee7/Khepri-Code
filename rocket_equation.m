function m_fuel = rocket_equation(delta_v,mf,Isp)
% PLEASE READ INPUTS BEFORE USE

%Inputs: Delta-v in m/s, final/dry mass of spacecraft (after fuel has 
%been burned) in kg, Isp in s is specific impulse to find amount
%of fuel
%Outputs: Mass of fuel in kg

g0 = 9.81;
m0 = mf*exp(delta_v/(Isp*g0));
m_fuel = m0-mf;
end
