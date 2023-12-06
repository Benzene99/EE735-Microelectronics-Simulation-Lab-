clc ;
clear vars ;
clear all;

%declaring given parameters
q = 1.6e-19;
Nd= 1e15;
Wp = 2e-6;
Wn = 2e-6;
epsilono = 8.85e-14;
T = 300 + 1.5*3; % roll number ends with 3
Vt = 0.02600 * (T/300);
Vbi = 0.6;
epsilonSi = 11.8;
ni = 1e10;

%calculating acceptor donor concentration

Na = (ni^2/Nd)*exp(Vbi/Vt);

%Depletion width calculation
Ld = sqrt(2*epsilonSi*epsilono*Vbi*q^(-1)*(Na^(-1)+Nd^(-1)));
xn = 10^4*(Na/(Na+Nd))*Ld;
xp = 10^4*(Nd/(Na+Nd))*Ld;

%Charge profile plot

dx = 0.0001; %step size assumed
P = -xp:dx:xn;
rho = zeros(1,length(P));


rho(1:7737) = (-q)*Na; 

rho(7738:13528) = (+q)*Nd; 

figure
plot( P * 10^(-4), rho)
hold on
xlabel('x(cm)');
ylabel('Volume Charge Density(C/cm^3)');
title('Volume Charge Density Profile');

%Electric field Profile
E= zeros(1,length(P));
%Trapezoidal Integration Method
for i=2:length(P)
    
    E(1,i) = E(1,i-1) + dx*0.5*(rho(1,i-1)+ rho(1,i));

end

% due to an arbitrary step size chosen (not in micrometers) we have to multiply 10e-4 to match the units and bring it in cm
E = 10^(-4) * (epsilonSi*epsilono)^(-1)*E;

figure
plot( P * 10^(-4), E)
hold on
xlabel('x(cm)');
ylabel('Electric Field(V/cm)');
title('Electric Field  Profile');

%Voltage Profile using trapezoidal rule
V= zeros(1,length(P));
for i=2:length(P)
    
    V(1,i) = V(1,i-1) + dx*0.5*(E(1,i-1)+ E(1,i));
   
end
% due to an arbitrary step size chosen (not in micrometers) we have to multiply 10e-4 to match the units and bring it in cm
V = -10^(-4)*V;

figure
plot( P * 10^(-4), V)
hold on
xlabel('x(cm)');
ylabel('Electric Potential(V)');
title('Electric Potential Profile');


