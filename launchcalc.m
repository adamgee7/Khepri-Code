function delta_v = launchcalc(n,u_eq,R) 

%Inputs
%n = number of serial stages
%u_eq = Isp * ge, if parallel, calculate an equivalent Isp by summing
       %thrusts of each engine and dividing by sum of their mass flow rates
%R = M0/Mb; %mass ratio, where M0 = Ml (payload mass) + Mp (propellant
     %mass) + Ms (structural mass), Mb = Ml + Ms

%Initialize
delta_v = 0;

%Calculate delta v
if n == 1 %single-stage
    delta_v = u_eq*log(R); %(m/s)
else %multi-stage
        for i=1:n
            delta_v = delta_v + u_eq(i)*log(R(i));
        end
end

end