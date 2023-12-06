clc ;
clear vars ;
clear all;

%declaring given parameters
q = 1.6e-19;
Nd= 1e15;
Wp = 2;
Wn = 2;
epsilono = 8.85e-14;
T = 300 + 1.5*3; % roll number ends with 3
Vt = 0.02600 * (T/300);
Vbi = 0.6;
epsilonSi = 11.8;
ni = 1e10;

%calculating acceptor donor concentration



%Depletion width calculation(in cm)
Ld = sqrt(6*epsilonSi*epsilono*Vbi*(q*ni)^(-1)*exp(-(Vbi/(2*Vt))));
xn = 10^4*0.5*Ld;
xp = 10^4*0.5*Ld;

a = 2*ni*exp(Vbi/(2*Vt))*Ld^(-1);
%Charge profile plot

dx = 0.0001; %step size assumed
X = -Wp:dx:Wn;
rho = zeros(1,length(X));

rho(1,round((-xp+Wp)/dx+1):round((xn+Wp)/dx+1)) = q*a*X(1,round((-xp+Wp)/dx+1):round((xn+Wp)/dx+1))*10^(-4);

figure
plot( X * 10^(-4), rho)
hold on
xlabel('x(cm)');
ylabel('Volume Charge Density(C/cm^3)');
title('Volume Charge Density Profile');

%Electric field Profile
E= zeros(1,length(X));
%Trapezoidal Integration Method
for i=2:length(X)
    
    E(1,i) = E(1,i-1) + dx*0.5*(rho(1,i-1)+ rho(1,i));

end

% due to an arbitrary step size chosen (not in micrometers) we have to multiply 10e-4 to match the units and bring it in cm
E = 10^(-4) * (epsilonSi*epsilono)^(-1)*E;

figure
plot( X * 10^(-4), E)
hold on
xlabel('x(cm)');
ylabel('Electric Field(V/cm)');
title('Electric Field  Profile');

%Voltage Profile using trapezoidal rule
V= zeros(1,length(X));
for i=2:length(X)
    
    V(1,i) = V(1,i-1) + dx*0.5*(E(1,i-1)+ E(1,i));
   
end
% due to an arbitrary step size chosen (not in micrometers) we have to multiply 10e-4 to match the units and bring it in cm
V = -10^(-4)*V;

figure
plot( X * 10^(-4), V)
hold on
xlabel('x(cm)');
ylabel('Electric Potential(V)');
title('Electric Potential Profile');
