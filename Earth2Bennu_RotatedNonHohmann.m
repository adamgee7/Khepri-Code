clc;
clear all;

%General Parameters
au = 1.496*10^11; %in km
m_E = 5.972 * 10^24;
m_S = 1.989 * 10^30;
R_S     = 6.957*10^5; % Sun radius (Km)
m_B = 7.329 * 10^10;
G = 6.6743 * 10^(-11);
mu = G*m_S;
theta_circle = 0:0.01:2*pi;

%Earth Orbit parameters:
e_E = 0.01671123; % Orbit eccentricity
i_E = 23.4393;
a_E = 1.000001018*au; %Semimajor axis (km)
r_E = a_E*(1-e_E^2)./(1+e_E*cos(theta_circle));
        
%Bennu orbit parameters
e_b = 0.20375; % Orbit eccentricity
i_b = 175; %degrees
a_b = 1.1264*au; %Semimajor axis (Km)

% Launch to DSM1
% Earth Orbit
e(1) = e_E;
r_A(1) = a_E*(1+e_E);
r_P(1) = a_E*(1-e_E);
h(1) = sqrt(2*mu*r_A(1)*r_P(1)/(r_A(1)+r_P(1)));
earth_rot = 0;
theta_start = -0.3706;
v(1) = mu/h(1)*sqrt(1+2*e(1)*cos(theta_start)+e(1)^2); %fastest at periapsis (nearing)
v_r(1) = mu/h(1)*e(1)*sin(theta_start);
v_n(1) = mu/h(1)*(1+e(1)*cos(theta_start));
gamma(1) = atan(v_r(1)/v_n(1));

% 1st End Orbit
theta_circle = 0:0.01:2*pi;
r_A(2) = r_A(1)*1.1;
r_P(2) = r_P(1)*0.85;
e(2) = (r_A(2)-r_P(2))/(r_A(2)+r_P(2));
h(2) = sqrt(2*mu*r_A(2)*r_P(2)/(r_A(2)+r_P(2)));
apse_angle(1)=80;
theta_rot(1) = earth_rot+apse_angle(1)*pi/180;
a = e(1)*h(2)^2-e(2)*h(1)^2*cos(theta_rot(1)-earth_rot);
b = -e(2)*h(1)^2*sin(theta_rot(1)-earth_rot);
c = h(1)^2-h(2)^2;
theta(1) = atan2(b,a)+acos(c/a*cos(atan2(b,a))); %wrt apse line of orbit 1
theta_offset(1) = theta(1)+earth_rot;
theta_end(1) = 2.0367;
v(2) = mu/h(2)*sqrt(1+2*e(2)*cos(theta_end(1)-theta_rot(1))+e(2)^2); %fastest at periapsis (nearing)
v_r(2) = mu/h(2)*e(2)*sin(theta_offset(1)-theta_rot(1));
v_n(2)= mu/h(2)*(1+e(2)*cos(theta_offset(1)-theta_rot(1)));
gamma(2) = atan(v_r(2)/v_n(2));
delta_v(1)=sqrt(v(1)^2+v(2)^2-2*v(1)*v(2)*cos(gamma(2)-gamma(1))); %Delta V from Launch (managed by Launch vehicle)

NonHohmann1(:,1) = theta_offset(1):0.01:theta_end(1);
a_NonHohmann(1) = (r_A(2)+r_P(2))/2;
NonHohmann1(:,2) = a_NonHohmann(1)*(1-e(2)^2)./(1+e(2)*cos(NonHohmann1(:,1)-theta_rot(1)));
NonHohmann1_full(:,1) = a_NonHohmann(1)*(1-e(2)^2)./(1+e(2)*cos(theta_circle-(theta_rot(1))));
v(3) = mu/h(2)*sqrt(1+2*e(2)*cos(theta_end(1)-theta_rot(1))+e(2)^2); %fastest at periapsis (nearing)
v_r(3) = mu/h(2)*e(2)*sin(theta_end(1)-theta_rot(1));
v_n(3) = mu/h(2)*(1+e(2)*cos(theta_end(1)-theta_rot(1)));
gamma(3) = atan(v_r(3)/v_n(3));
 
%DSM1 to Gravity Assist
r_A(3) = r_A(2)*1;
r_P(3) = r_A(2)*0.76;
e(3) = (r_A(3)-r_P(3))/(r_A(3)+r_P(3));
h(3) = sqrt(2*mu*r_A(3)*r_P(3)/(r_A(3)+r_P(3)));
apse_angle(2)=95;
theta_rot(2) =theta_rot(1)+(apse_angle(2)-apse_angle(1))*pi/180; %wrt apse line of previous orbit
a = e(2)*h(3)^2-e(3)*h(2)^2*cos(theta_rot(2)-theta_rot(1));
b = -e(3)*h(2)^2*sin(theta_rot(2)-theta_rot(1));
c = h(2)^2-h(3)^2;
theta(2) = atan2(b,a)+acos(c/a*cos(atan2(b,a))); %wrt apse line of orbit 1
theta_offset(2) = theta(2)+theta_rot(1);
theta_end(2) = 2*pi-0.0271;
v(4) = mu/h(3)*sqrt(1+2*e(3)*cos(theta_offset(2)-theta_rot(2))+e(3)^2);
v_r(4) = mu/h(3)*e(3)*sin(theta_offset(2)-theta_rot(2));
v_n(4) = mu/h(3)*(1+e(3)*cos(theta_offset(2)-theta_rot(2)));
gamma(4) = atan(v_r(4)/v_n(4));
delta_v(2)=sqrt(v(3)^2+v(4)^2-2*v(3)*v(4)*cos(gamma(4)-gamma(3))); %DSM1 DeltaV

NonHohmann2(:,1) = theta_offset(2):0.01:theta_end(2);
a_NonHohmann(2) = (r_A(3)+r_P(3))/2;
NonHohmann2(:,2) = a_NonHohmann(2)*(1-e(3)^2)./(1+e(3)*cos(NonHohmann2(:,1)-theta_rot(2)));
NonHohmann2_full(:,1) = a_NonHohmann(2)*(1-e(3)^2)./(1+e(3)*cos(theta_circle-(theta_rot(2))));
v(5) = mu/h(3)*sqrt(1+2*e(3)*cos(theta_end(2)-theta_rot(2))+e(3)^2);
v_r(5) = mu/h(3)*e(3)*sin(theta_end(2)-theta_rot(2));
v_n(5) = mu/h(3)*(1+e(3)*cos(theta_end(2)-theta_rot(2)));
gamma(5) = atan(v_r(5)/v_n(5));

%Gravity Assist to DSM2
r_A(4) = r_A(3)*1.18;
r_P(4) = r_A(3)*0.79;
e(4) = (r_A(4)-r_P(4))/(r_A(4)+r_P(4));
h(4) = sqrt(2*mu*r_A(4)*r_P(4)/(r_A(4)+r_P(4)));
apse_angle(3)=65;
theta_rot(3) =theta_rot(2)+(apse_angle(3)-apse_angle(2))*pi/180; %wrt apse line of previous orbit
a = e(3)*h(4)^2-e(4)*h(3)^2*cos(theta_rot(3)-theta_rot(2));
b = -e(4)*h(3)^2*sin(theta_rot(3)-theta_rot(2));
c = h(3)^2-h(4)^2;
theta(3) = atan2(b,a)+acos(c/a*cos(atan2(b,a))); %wrt apse line of orbit 1
theta_offset(3) = theta(3)+theta_rot(2)-2*pi;
theta_end(3) = 2*pi-2.1334;
v(6) = mu/h(4)*sqrt(1+2*e(4)*cos(theta_offset(3)-theta_rot(3))+e(4)^2);
v_r(6) = mu/h(4)*e(4)*sin(theta_offset(3)-theta_rot(3));
v_n(6) = mu/h(4)*(1+e(4)*cos(theta_offset(3)-theta_rot(3)));
gamma(6) = atan(v_r(6)/v_n(6));
delta_v(3)=sqrt(v(5)^2+v(6)^2-2*v(5)*v(6)*cos(gamma(6)-gamma(5))); %Gravity Assist DeltaV (done by Earth)

NonHohmann3(:,1) = theta_offset(3):0.01:theta_end(3);
a_NonHohmann(3) = (r_A(4)+r_P(4))/2;
NonHohmann3(:,2) = a_NonHohmann(3)*(1-e(4)^2)./(1+e(4)*cos(NonHohmann3(:,1)-theta_rot(3)));
NonHohmann3_full(:,1) = a_NonHohmann(3)*(1-e(4)^2)./(1+e(4)*cos(theta_circle-(theta_rot(3))));
v(7) = mu/h(4)*sqrt(1+2*e(4)*cos(theta_end(3)-theta_rot(3))+e(4)^2);
v_r(7) = mu/h(4)*e(4)*sin(theta_end(3)-theta_rot(3));
v_n(7) = mu/h(4)*(1+e(4)*cos(theta_end(3)-theta_rot(3)));
gamma(7) = atan(v_r(7)/v_n(7));

%DSM2 to Rendezvous
r_A(5) = r_A(4)*1;
r_P(5) = r_A(4)*0.7;
e(5) = (r_A(5)-r_P(5))/(r_A(5)+r_P(5));
h(5) = sqrt(2*mu*r_A(5)*r_P(5)/(r_A(5)+r_P(5)));
apse_angle(4)=50;
theta_rot(4) =theta_rot(3)+(apse_angle(4)-apse_angle(3))*pi/180; %wrt apse line of previous orbit
a = e(4)*h(5)^2-e(5)*h(4)^2*cos(theta_rot(4)-theta_rot(3));
b = -e(5)*h(4)^2*sin(theta_rot(4)-theta_rot(3));
c = h(4)^2-h(5)^2;
theta(4) = atan2(b,a)+acos(c/a*cos(atan2(b,a))); %wrt apse line of orbit 1
theta_offset(4) = theta(4)+theta_rot(3)-2*pi;
theta_end(4) = 0.0989;
v(8) = mu/h(5)*sqrt(1+2*e(5)*cos(theta_offset(4)-theta_rot(4))+e(5)^2);
v_r(8) = mu/h(5)*e(5)*sin(theta_offset(4)-theta_rot(4));
v_n(8) = mu/h(5)*(1+e(5)*cos(theta_offset(4)-theta_rot(4)));
gamma(8) = atan(v_r(8)/v_n(8));
delta_v(4)=sqrt(v(7)^2+v(8)^2-2*v(7)*v(8)*cos(gamma(8)-gamma(7))); %DSM2 DeltaV

NonHohmann4(:,1) = theta_offset(4):0.01:theta_end(4);
a_NonHohmann(4) = (r_A(5)+r_P(5))/2;
NonHohmann4(:,2) = a_NonHohmann(4)*(1-e(5)^2)./(1+e(5)*cos(NonHohmann4(:,1)-theta_rot(4)));
NonHohmann4_full(:,1) = a_NonHohmann(4)*(1-e(5)^2)./(1+e(5)*cos(theta_circle-(theta_rot(4))));
v(9) = mu/h(5)*sqrt(1+2*e(5)*cos(theta_end(4)-theta_rot(4))+e(5)^2);
v_r(9) = mu/h(5)*e(5)*sin(theta_end(4)-theta_rot(4));
v_n(9) = mu/h(5)*(1+e(5)*cos(theta_end(4)-theta_rot(4)));
gamma(9) = atan(v_r(9)/v_n(9));

% Bennu Orbit
e(6) = e_b;
r_A(6) = a_b*(1+e_b);
r_P(6) = a_b*(1-e_b);
h(6) = sqrt(2*mu*r_A(6)*r_P(6)/(r_A(6)+r_P(6)));
bennu_rot = 60;
r_b = a_b*(1-e_b^2)./(1+e_b*cos(theta_circle-bennu_rot*pi/180)); %Offset to match OSIRIS REx diagram
theta_rot(5) =theta_rot(4)+(bennu_rot-apse_angle(4))*pi/180; %apse line at 40 deg
a = e(5)*h(6)^2-e(6)*h(5)^2*cos(theta_rot(5)-theta_rot(4));
b = -e(6)*h(5)^2*sin(theta_rot(5)-theta_rot(4));
c = h(5)^2-h(6)^2;
theta(5) = atan2(b,a)+acos(c/a*cos(atan2(b,a))); %wrt apse line of orbit 1
theta_offset(5) = theta(5)+theta_rot(4)
v(10) = mu/h(6)*sqrt(1+2*e(6)*cos(theta_offset(5)-theta_rot(5))+e(6)^2);
v_r(10) = mu/h(6)*e(6)*sin(theta_offset(5)-theta_rot(5));
v_n(10) = mu/h(6)*(1+e(6)*cos(theta_offset(5)-theta_rot(5)));
gamma(10) = atan(v_r(10)/v_n(10));
delta_v(5)=sqrt(v(9)^2+v(10)^2-2*v(9)*v(10)*cos(gamma(10)-gamma(9))) %Rendezvous DeltaV
total_spacecraft_dv = delta_v(2)+delta_v(4)+delta_v(5) %Total impulses of Spacecraft

% Plotting
hpp1 = polarplot(theta_circle,r_E,'g');
hold on
%polarplot([0 0+180; 0 0+180]*pi/180, [0 0; 1 1]*a_E*1.5,'--r')

hpp2 = polarplot(NonHohmann1(:,1),NonHohmann1(:,2),'r');
%hpp1 = polarplot(theta_circle,NonHohmann1_full(:,1),'k');
%polarplot([theta_rot(1)*180/pi theta_rot(1)*180/pi+180; theta_rot(1)*180/pi theta_rot(1)*180/pi+180]*pi/180, [0 0; 1 1]*a_E*1.5,'--k')
polarplot(theta_offset(1),NonHohmann1(1,2),'or')

hpp3 = polarplot(NonHohmann2(:,1),NonHohmann2(:,2),'r');
%polarplot(theta_circle,NonHohmann2_full(:,1),'g')
%polarplot([theta_rot(2)*180/pi theta_rot(2)*180/pi+180; theta_rot(2)*180/pi theta_rot(2)*180/pi+180]*pi/180, [0 0; 1 1]*a_E*1.5,'--g')
polarplot(theta_offset(2),NonHohmann2(1,2),'or')

hpp4 = polarplot(NonHohmann3(:,1),NonHohmann3(:,2),'r');
%polarplot(theta_circle,NonHohmann3_full(:,1),'m')
%polarplot([theta_rot(3)*180/pi theta_rot(3)*180/pi+180; theta_rot(3)*180/pi theta_rot(3)*180/pi+180]*pi/180, [0 0; 1 1]*a_E*1.5,'--m')
polarplot(theta_offset(3),NonHohmann3(1,2),'or')

hpp5 = polarplot(NonHohmann4(:,1),NonHohmann4(:,2),'r');
%polarplot(theta_circle,NonHohmann4_full(:,1),'c')
%polarplot([theta_rot(4)*180/pi theta_rot(4)*180/pi+180; theta_rot(4)*180/pi theta_rot(4)*180/pi+180]*pi/180, [0 0; 1 1]*a_E*1.5,'--c')
polarplot(theta_offset(4),NonHohmann4(1,2),'or')

hpp6 = polarplot(theta_circle,r_b,'b');
hold on
%polarplot([bennu_rot bennu_rot+180; bennu_rot bennu_rot+180]*pi/180, [0 0; 1 1]*a_E*1.5,'--b')
polarplot(theta_offset(5),NonHohmann4(end,2),'or')
legend([hpp1(1),hpp2(1),hpp3(1),hpp4(1),hpp5(1),hpp6(1)],'Earth','Launch to DSM1','DSM1 to Gravity Assist','Gravity Assist to DSM2','DSM2 to Rendezvous','Bennu')
