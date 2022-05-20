function [time,h_max,arc_length,lateral_length] = ejection_calc(alpha,v)
% PLEASE READ INPUTS BEFORE USE

%Inputs: alpha is the angle relative to the surface that the particle is 
% ejected at (in degrees); v is the ejection velocity (initial velocity) in
% m/s, arc_length and lateral_length are the respective lengths in m
%Outputs: time in seconds, maximum height of the projectile

% Constants and Parameters
r_Bennu = 245.03; %m from Wikipedia
G = 6.67*10^(-11);
M_Bennu = 7.329 * 10^(10); %kg from Wikipedia
mu = G*M_Bennu;

% Particle Ejection Calculator
alpha = alpha*pi/180;
h = r_Bennu*v*cos(alpha);
theta = atan2(sin(alpha)*h^2,cos(alpha)*(h^2-mu*r_Bennu));
e = sin(alpha)/sin(theta-alpha);

if v > sqrt(2*mu/r_Bennu)
    a = h^2/mu*1/(e^2-1);
    h_max = NaN;
    time = NaN;
    arc_length = NaN;
    lateral_length = NaN;
    
    %Plotting
    theta_circle = 0:0.01:2*pi;
    rb_circle = r_Bennu*ones(length(theta_circle));
    theta_projectile = theta:0.0001:theta+(pi-theta)/2;
    theta_other = 0:0.001:theta;
    r_projectile = a*(1-e^2)./(1+e*cos(theta_projectile));
    r_other = a*(1-e^2)./(1+e*cos(theta_other));
    hpp1 = polarplot(theta_circle,rb_circle,'r');
    rlim([0 700])
    hold on
    hpp2 = polarplot(theta_projectile,r_projectile,'k');
    hpp3 = polarplot(theta_other,r_other,'--k');
    legend([hpp1(1),hpp2(1)],'Bennu','Projectile');
elseif alpha == 90
    sprintf("Retry with a smaller angle. This particle cannot be ejected.")
    time = NaN;
    h_max = NaN;
    arc_length = NaN;
    lateral_length = NaN;
elseif alpha == 0
    sprintf("Retry with a larger angle. This particle cannot be ejected.")
    time = NaN;
    h_max = NaN;
    arc_length = NaN;
    lateral_length = NaN;
else
    a = h^2/mu*1/(1-e^2);
    r_a = a*(1+e);
    h_max = r_a-r_Bennu;
    E = 2*atan(sqrt((1-e)/(e+1))*tan(theta/2*pi/180));
    M_e = E-e*sin(E);
    time = h^3/mu^2/sqrt(1-e^2)^3*M_e*2; %times 2 for projectile, seems awfully short.... ? or else backwards?
    
    fun = @(x) a*sqrt(1-e^2.*cos(x).^2);
    arc_length = integral(fun,theta,2*pi-theta);
    lateral_length = r_Bennu*sin(pi-theta)*2;
    
    %Plotting
    theta_circle = 0:0.01:2*pi;
    rb_circle = r_Bennu*ones(length(theta_circle));
    theta_projectile = theta:0.0001:(2*pi-theta);
    theta_other = 0:0.001:theta+(2*pi-theta):0.01:2*pi;
    r_projectile = a*(1-e^2)./(1+e*cos(theta_projectile));
    r_other = a*(1-e^2)./(1+e*cos(theta_other));
    
    hpp1 = polarplot(theta_circle,rb_circle,'r');
    hold on
    hpp2 = polarplot(theta_projectile,r_projectile,'k');
    hpp3 = polarplot(theta_other,r_other,'--k');
    legend([hpp1(1),hpp2(1)],'Bennu','Projectile')
end
end
