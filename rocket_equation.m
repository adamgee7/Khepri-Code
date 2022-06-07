function m_fuel = rocket_equation(delta_v,mf,Isp)
g0 = 9.81;
m0 = mf*exp(delta_v/(Isp*g0));
m_fuel = m0-mf;
end
