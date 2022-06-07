function [delta_v, m_fuel, delta_energy, time, orbit, orbital_period] = bennu2orbit(height,mode,mf,theta2,plot_cond,Isp) 
% PLEASE READ INPUTS BEFORE USE

%Inputs: height above surface; mode = hohmann, non-hohmann/hyperbolic, 
%escape or lob; theta2 = angle of departure for non-hohmann, plot_cond = 1
%to plot or anything else otherwise, Isp in s is specific impulse to find amount
%of fuel
%Outputs: Delta-v in m/s, specific change in energy in J, transit time in 
%hours, orbit type

% Notes: This calculator has 4 options for orbital maneuvers around Bennu.
% Two impulse maneuvers like Hohmann Transfers are the most energy
% efficient, compared to two impulse Non-Hohmann Transfers. One impulse
% maneuvers set objects into an elliptical orbit from the surface, or on an
% escape trajectory to leave Bennu.

% Constants and Parameters
r_Bennu = 245.03; %m from Wikipedia
G = 6.67*10^(-11);
M_Bennu = 7.329 * 10^(10); %kg from Wikipedia
mu = G*M_Bennu;
au = 1.496*10^11; %in m
a_Bennu = 1.1264*au;
M_Sun = 1.989 * 10^30;

% Orbital Maneuver Calculations
r1 = r_Bennu;
r2 = r_Bennu+height;
orbital_period = 2*pi/sqrt(mu)*((r2+r1)/2)^(3/2)/3600;
r_H = a_Bennu*(M_Bennu/(3*M_Sun))^(1/3); %Hill sphere radius around Bennu in m
if r2 < r_H
    % Two impulse maneuvers
    if mode == "Hohmann"
        v1 = sqrt(mu/r1);
        v_transfer_1 = sqrt(mu)*(sqrt(2/r1-2/(r1+r2)));
        delta_v_1 = v_transfer_1 - v1;
        delta_energy_1 = v_transfer_1^2/2-v1^2/2;
    
        v_transfer_2 = sqrt(mu)*(sqrt(2/r2-2/(r1+r2)));
        v2 =  sqrt(mu/r2);
        delta_v_2 = v2 - v_transfer_2;
        delta_energy_2 = v2^2/2-v_transfer_2^2/2;
    
        delta_v = delta_v_1 + delta_v_2;
        delta_energy = delta_energy_1 + delta_energy_2;
        time = pi/sqrt(mu)*((r2+r1)/2)^(3/2)/3600; % Time in hours
        orbit = "elliptic";

        %Plotting
        if plot_cond == 1
            theta_circle = 0:0.01:2*pi;
            r1_circle = r1*ones(length(theta_circle));
            r2_circle = r2*ones(length(theta_circle));
            theta_Hohmann = 0:0.01:pi;
            e_Hohmann = (r2-r1)/(r1+r2);
            a_Hohmann = r2/(1+e_Hohmann);
            r_Hohmann = a_Hohmann*(1-e_Hohmann^2)./(1+e_Hohmann*cos(theta_Hohmann));
    
            hpp1 = polarplot(theta_circle,r1_circle,'r');
            hold on
            hpp2 = polarplot(theta_circle,r2_circle,'b');
            hpp3 = polarplot(theta_Hohmann,r_Hohmann,'k');
            %legend('Bennu','Final Orbit')
            hold off
            legend([hpp1(1),hpp2(1),hpp3(1)],'Bennu','Orbit','Transfer')
        end

    end
    if mode == "Non-Hohmann" 
        if theta2 <= 180 % This angle depends on the radius.... - make a condition where if a NAN (ie not possible) then print saying to increase the angle
            %On Bennu
            h1 = sqrt(mu*r1);
            e1 = 0; %circular orbit on Bennu surface
            theta1 = 0;
            v1 = mu/h1*(1+e1*cos(theta1*pi/180));
            fpa_1 = 0;
            
            %Transfer Orbit
            e_transfer = (r2-r1)/(r1*cos(theta1*pi/180)-r2*cos(theta2*pi/180));
            h_transfer = sqrt(mu*r1*r2)*sqrt((cos(theta2*pi/180)-cos(theta1*pi/180))/(r2*cos(theta2*pi/180)-r1*cos(theta1*pi/180)));
            v_r_transfer_1 = mu/h_transfer*e_transfer*sin(theta1*pi/180);
            v_n_transfer_1 = mu/h_transfer*(1+e_transfer*cos(theta1*pi/180));
            v_transfer_1 = sqrt(v_r_transfer_1^2+v_n_transfer_1^2);
            fpa_transfer_1 = 180/pi*atan(v_r_transfer_1/v_n_transfer_1); %flight path angle of departure
            delta_fpa_1 = fpa_transfer_1-fpa_1;
            delta_v_1 = sqrt(v1^2+v_transfer_1^2-2*v1*v_transfer_1*cos(delta_fpa_1*pi/180));
            delta_energy_1 = v_transfer_1^2/2-v1^2/2;
            
            v_r_transfer_2 = mu/h_transfer*e_transfer*sin(theta2*pi/180);
            v_n_transfer_2 = mu/h_transfer*(1+e_transfer*cos(theta2*pi/180));
            v_transfer_2 = sqrt(v_r_transfer_2^2+v_n_transfer_2^2);
            fpa_transfer_2 = 180/pi*atan(v_r_transfer_2/v_n_transfer_2); 
        
            %Final Orbit
            h2 = sqrt(mu*r2); %circular orbit at height above Bennu surface
            e2 = 0;
            v2 = mu/h2*(1+e2*cos(theta2*pi/180));
            fpa_2 = 0;
            delta_fpa_2 = fpa_2-fpa_transfer_2;
            delta_v_2 = sqrt(v_transfer_2^2+v2^2-2*v2*v_transfer_2*cos(delta_fpa_2*pi/180));
            delta_energy_2 = v2^2/2-v_transfer_2^2/2;
        
            %Overall
            delta_v = delta_v_1 + delta_v_2;
            delta_energy = delta_energy_1 + delta_energy_2;
        
            %Time
            if e_transfer < 1 && e_transfer > 0
                E = 2*atan(sqrt((1-e_transfer)/(e_transfer+1))*tan(theta2/2*pi/180));
                M_e = E-e_transfer*sin(E);
                time = h_transfer^3/mu^2/sqrt(1-e_transfer^2)^3*M_e/3600;
                orbit = "elliptic";
                
                if plot_cond == 1
                    %Plotting
                    theta_circle = 0:0.01:2*pi;
                    r1_circle = r1*ones(length(theta_circle));
                    r2_circle = r2*ones(length(theta_circle));
                    theta_NonHohmann = 0:0.01:theta2*pi/180;
                    a_NonHohmann = r1/(1-e_transfer);
                    r_NonHohmann = a_NonHohmann*(1-e_transfer^2)./(1+e_transfer*cos(theta_NonHohmann));
            
                    hpp1 = polarplot(theta_circle,r1_circle,'r');
                    hold on
                    hpp2 = polarplot(theta_circle,r2_circle,'b');
                    hpp3 = polarplot(theta_NonHohmann,r_NonHohmann,'k');
                    legend([hpp1(1),hpp2(1),hpp3(1)],'Bennu','Orbit','Transfer')
                end

            elseif e_transfer > 1
                E_bar = 2*atanh(sqrt((e_transfer-1)/(e_transfer+1))*tan(theta2/2*pi/180));
                M_h = e_transfer*sinh(E_bar)-E_bar;
                time = h_transfer^3/mu^2/sqrt(e_transfer^2-1)^3*M_h/3600;
                orbit = "hyperbolic";
                
                if plot_cond == 1
                    %Plotting
                    theta_circle = 0:0.01:2*pi;
                    r1_circle = r1*ones(length(theta_circle));
                    r2_circle = r2*ones(length(theta_circle));
                    theta_NonHohmann = 0:0.01:theta2*pi/180;
                    a_NonHohmann = r1/(e_transfer-1);
                    r_NonHohmann = a_NonHohmann*(e_transfer^2-1)./(1+e_transfer*cos(theta_NonHohmann));
            
                    hpp1 = polarplot(theta_circle,r1_circle,'r');
                    hold on
                    hpp2 = polarplot(theta_circle,r2_circle,'b');
                    hpp3 = polarplot(theta_NonHohmann,r_NonHohmann,'k');
                    legend([hpp1(1),hpp2(1),hpp3(1)],'Bennu','Orbit','Transfer')
                end

            elseif e_transfer < 0
                sprintf("Retry with a larger angle but less than 180 degrees. This angle deparature is impossible.")
                delta_v = NaN;
                delta_energy = NaN;
                time = NaN;
                orbit = "None";
            end
        else
            sprintf("Retry with an angle less than 180 degrees.")
                delta_v = NaN;
                delta_energy = NaN;
                time = NaN;
                orbit = "None";
        end
    end
    
    % One impulse Maneuvers
    if mode == "Escape"
        v1 = sqrt(mu/r1);
        v2 = sqrt(2*mu/r_Bennu);
        delta_v = v2-v1;
        time = NaN;
        delta_energy = v2^2/2-v1^2/2;
        orbit = "hyperbolic";
        
        if plot_cond == 1
            %Plotting
            e_escape = 1.00001; %This is a function of the ejection speed
            a_escape = r1/(e_escape-1);
            v_hyp = sqrt(mu*(1+2*e_escape+e_escape^2)/(a_escape*(e_escape^2-1)));
            theta_circle = 0:0.01:2*pi;
            r1_circle = r1*ones(length(theta_circle));
            theta_escape = 0:0.01:100*pi/180;
            r_escape = a_escape*(e_escape^2-1)./(1+e_escape*cos(theta_escape));
    
            hpp1 = polarplot(theta_circle,r1_circle,'r');
            hold on
            hpp2 = polarplot(theta_escape,r_escape,'k');
            legend([hpp1(1),hpp2(1)],'Bennu','Escape Trajectory')
        end

    end
    if mode == "Lob"
        v1 = sqrt(mu/r1);
        e = (r2-r1)/(r1+r2);
        v2 = sqrt(mu/r1*(1+e));
        delta_v = v2-v1;
        delta_energy = v2^2/2-v1^2/2;
        time = pi/sqrt(mu)*((r2+r1)/2)^(3/2)/3600;
        orbit = "elliptic";
        
        if plot_cond == 1
            %Plotting
            theta_circle = 0:0.01:2*pi;
            r1_circle = r1*ones(length(theta_circle));
            a_Lob = r1/(1-e);
            r_Lob = a_Lob*(1-e^2)./(1+e*cos(theta_circle));
    
            hpp1 = polarplot(theta_circle,r1_circle,'r');
            hold on
            hpp2 = polarplot(theta_circle,r_Lob,'k');
            legend([hpp1(1),hpp2(1)],'Bennu','Lob Orbit')
        end
    end
m_fuel = rocket_equation(delta_v,mf,Isp);
else
    sprintf("Outside of Bennu's gravitational sphere of influence (Hill Sphere). Try a smaller height.")
    delta_v = NaN;
    delta_energy = NaN;
    time = NaN;
    orbit = "None";
end
end