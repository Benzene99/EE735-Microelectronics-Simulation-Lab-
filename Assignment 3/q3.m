clc,
clear vars
clear all;

% Defining Constants
L = 4e-4;% Length of n+-n junction diode with symmetric width on both sides taken
X_num = 300;% number of grid points along x axis considered
dx = L/X_num; %grid spacing calculated
Vbi = 0.06+(3*3/300); % Roll no is 23M1133
ni = 1.34e10;%intrinsic carrier conc to match Vbi value within 5% tolerance
T=300;%Temperature
Vth = 0.026;%Thermal voltage
Nd2 = 1.18e14;% selected to match Vbi value within 5% tolerance
q = 1.6e-19;
epslon = 11.8*8.85e-14 ;
nL= 2e-4;
Nd1 = Nd2*exp(Vbi/Vth); % Nd1 chosen to match Vbi value within 5% tolerance

%Initialising Voltage,charge,carrier concentration matrices 
V = zeros(1,X_num);
V1 = zeros(1,X_num);
Ef = zeros(1,X_num);
rho = zeros(1,X_num);
p = zeros(1,X_num);
n = zeros(1,X_num);
diagBV = zeros(X_num,X_num);

% Defining A matrix to use it in LU
A = zeros(X_num);        
for i=1:X_num
    for j=1:X_num
        if i==1&&j==1
            A(i,j)=1;
        elseif i==X_num&&j==X_num
            A(i,j)=1;
        elseif i==j
            A(i,j)=-2;
        elseif (i==1&&j==2)||(i==X_num&&j==X_num-1)
            A(i,j)=-1;
        elseif i==j-1
            A(i,j)=1;
        elseif i==j+1
            A(i,j)=1;
        end
    end
end


for z = 1:500 %no of iterations

    V(1,X_num) = Vth*log(Nd2/ni) ;
    V(1,1) = Vth*log(Nd1/ni);
    p = ni*exp(-V/Vth);
    n = ni*exp(V/Vth);
    rho(1,1) = 0 ;
    rho(1,2:X_num/2) = p(1,2:X_num/2)-n(1,2:X_num/2)+Nd1 ;
    rho(1,X_num/2+1:X_num-1)= p(1,X_num/2+1:X_num-1)-n(1,X_num/2+1:X_num-1)+Nd2;
    rho(1,X_num) = 0 ; 
    B = ((-q*dx^2)/epslon)*rho';
    F = A*V'-B;
    diagBV = diag((((-q*dx^2)/epslon)*(-ni/Vth)*(exp(-V/Vth)+exp(V/Vth))));

    J = A-diagBV;
    %LU Decompistion  

    a = J;    
    x = V1';
    b = F;   %charge Density RHS of Poisson equation

    % Diagonalization
    for k = 1 : X_num-1
    for i = k+1 : X_num
    a(i,k) = a(i,k)/a(k,k);
    for j = k+1:X_num
    a(i,j) = a(i,j)-a(i,k)*a(k,j);
    end
    end
    end
    % Forward substitution to solve Ld=b
    x(1)=b(1) ;
    for i = 2 : X_num
    s=0;
    for j = 1 : i-1
    s = s + a(i,j)*x(j);
    end
    x(i) = b(i)-s;
    end
    % Backward substitution to solve Ux=d
    x(X_num) = x(X_num) / a(X_num, X_num) ;
    for i = X_num-1 : -1 : 1
    s=0 ;
    for j = i+1:X_num
    s = s + a( i, j) * x(j) ;
    end
    x(i) = (x(i) - s)/a(i,i) ;
    end
    V1 = x';
    V = V - V1;
end


for x = 1:X_num
    Ef(1,x) = (((X_num-1)*dx)/(X_num-2))*(q/epslon)*(sum(rho(1,1:x)));  % using Trapezoidal integral method to Calulate the Electric Field
end


figure
plot(1:X_num,V,'DisplayName','Actual Voltage Profile',"LineWidth",2);
xlabel('Distance(along X axis)');
ylabel('Electric Potential') ;
title('Electric Potential Vs Distance') ;
grid ;
legend;


figure
plot(1:X_num,rho,'DisplayName','Actual Charge Profile',"LineWidth",2);
xlabel('Distance(along X axis)');
ylabel('ChargeDensity') ;
title('ChargeDensity Vs Distance') ;
grid;
legend ;


figure
semilogy(1:X_num,p,'DisplayName','Holes',"LineWidth",2);
xlabel('Distance(along X axis)');
ylabel('p or n') ;
title('Carrier Concentration Vs Distance') ;
grid ;
hold on;
semilogy(1:X_num,n,'DisplayName','Electrons',"LineWidth",2);
legend;
hold off

figure
plot(1:X_num,Ef,'DisplayName','Electric Field',"LineWidth",2);
xlabel('Distance(along X axis)');
ylabel('Electric Field') ;
title('Electric Field Vs Distance') ;
grid ;
legend;

Ef = zeros(1,X_num);
Ei = Ef - V ;
Eg= 1.12 ; %Silicon Band gap
Ec = Eg/2 + Ei ;
Ev = Ei-Eg/2 ;
figure
plot(1:X_num,Ec,'DisplayName','Ec',"LineWidth",2);
hold on
plot(1:X_num,Ei,'DisplayName','Ei',"LineWidth",2);
plot(1:X_num,Ef,'DisplayName','Ef',"LineWidth",2);
plot(1:X_num,Ev,'DisplayName','Ev',"LineWidth",2);
title('Energy Band Diagram') ;
legend;
hold off

fprintf('The lenth of N+ region is %d centimeters\n',L/2);
fprintf('The lenth of N region is from %d centimeters starting from %d centimeters\n',L/2,nL);
fprintf('The doping of N+ region is %d cm-3\n',Nd1);
fprintf('The doping of N region is %d cm-3\n',Nd2);
