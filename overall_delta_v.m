function [m_fuel_e2b, t_e2b, m_fuel_prox, t_prox] = overall_delta_v(m_pod,m_mothership,num_pods,m_plant,num_plant,orb_height,num_trips)

% Earth 2 Bennu
m_dry_e2b = m_pod*num_pods + m_plant*num_plant + m_mothership; %not including any fuel for mothership
m_OREx_dry = 880;
m_OREx_wet = 2110;
Isp_OREx(1) = 225;
Isp_OREx(2) = 236;
delta_v_OREx = log(m_OREx_wet/m_OREx_dry)*mean(Isp_OREx)*9.81;
SF = round(delta_v_OREx/1400,1);
delta_v_e2b = delta_v_OREx/SF-10*SF; %amount for half of O-REx getting rid of safety factor
m_fuel_e2b = rocket_equation(delta_v_e2b,m_dry_e2b);
t_e2b = 712; %days

% Proximity Ops at Bennu
theta2(1) = 180;
[delta_v(1), ~, ~, ~, ~, ~] = bennu2orbit(orb_height,"Non-Hohmann",0,theta2(1),0);
theta2(2) = 170;
i=2;
while 1
    [delta_v(i), ~, ~, ~, ~, ~] = bennu2orbit(orb_height,"Non-Hohmann",0,theta2(i),0);
    if ~isnan(delta_v(i))
        if isnan(delta_v(i-1))
            theta_end = theta2(i);
            break
        end
        theta2(i+1) = theta2(i)-10;
    else
        theta2(i+1) = theta2(i)+1;
    end
    i=i+1;
end
theta_array = linspace(180,theta_end,11);

for i=1:11
    [delta_v_b2o(i), ~, ~, t_b2o(i), ~, ~] = bennu2orbit(orb_height,"Non-Hohmann",0,theta_array(i),0);
    m_fuel_b2o_trip(i) = rocket_equation(delta_v_b2o(i)*2,m_pod);
end
m_fuel_prox = num_trips*m_fuel_b2o_trip;
t_prox = num_trips*t_b2o/24;
end