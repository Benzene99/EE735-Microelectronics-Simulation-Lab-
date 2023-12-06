% Q1.a- Abrupt PN Juntion
clc;
clear vars;
clear all;

%Defining all constants required for solving Poisson Equation
epsilono = 8.84e-12;
epsilonsi = 11.8* epsilono;
E_g = 1.12; 
q = 1.6e-19;
ni = 1e10 * 1e6;
Vth = 0.026; 
L = 2e-6; %lenth of PN Junction
dx = 1e-9; %step size for computation of AV=B optimised to get best results

dist = -L/2:dx:L/2;

% Calculating total grid points
x_num = round(L/dx)+1;

% Doping concentration according to Roll no

Na = 4e15 * 1e6; % Roll no is 23M1133 so a=3
Nd = 4e15 * 1e6; % Roll no is 23M1133 so b=3

% Boundary Condition (Potential Values)
V1 = -Vth * log(Na/ni);
%V1 = 0; 
V2 = Vth * log(Nd/ni);

% Initialising Doping Profile
Doping = zeros(x_num,1);
Doping(1:round(x_num/2)) = -Na;
Doping(round(x_num/2)+1:x_num) = Nd;

% Generating A matrix to solve Poisson Equation
A = zeros(x_num, x_num);
A(1,1) = 1;
A(x_num,x_num) = 1;
for i = 2:x_num-1
    A(i,i) = -2;
    A(i,i-1) = 1;
    A(i,i+1) = 1;
end

%Generating B matrix 
B = zeros(x_num,1);
%Applying boundary conditions given in the problem
B(1,1) = V1;
B(x_num,1) = V2;

% Initialize Voltage Matrix with all elements to zero except at the extreme boundary
% points 
V = zeros(x_num,1);
V(1,1) = V1;
V(x_num,1) = V2;

% Newton Raphson Iterations to solve V Matrix

err = 3.14e-13; % Error Tolerance for Newton Raphson Method

z = 0; %iteration counter
while(true)
    z = z+1;

    p = ni .* exp(-V./Vth);
    n = ni .* exp(V./Vth);
    
    for i = 2:round(x_num/2)
        B(i,1) = -q*(Doping(i) + p(i))*(dx^2 / epsilonsi);
    end
    for i = round(x_num/2)+1:x_num-1
        B(i,1) = -q*(Doping(i) - n(i))*(dx^2 / epsilonsi);
    end

    % Finding the jacobian matrix
    B_V = zeros(x_num,1);
    for i = 2:round(x_num/2)
        B_V(i,1) = q*p(i)*(dx^2)/(epsilonsi*Vth);
    end
    for i = round(x_num/2)+1:x_num-1
        B_V(i,1) = q*n(i)*(dx^2)/(epsilonsi*Vth);
    end

    Jac = A - diag(B_V);

    % Calculating f Matrix
    f = A*V - B;
    V_next = V - Jac\f;

    sim_error = max(V_next)-max(V); %quantifying simulated error

    if(abs(sim_error) < err)
        break;
    end
    V = V_next;
end

% Electric Field using gradient function
E = -gradient(V,dx);

% Plots
figure(1);
plot(dist,V,"LineWidth",2);
xlabel("Distance (along X axis)");
ylabel("Potential (Volts)");
title("Potential along X axis (Abrupt PN Junction)");

figure(2);
plot(dist,E,"LineWidth",2);
xlabel("Distance (along X axis)");
ylabel("Electric Field (Volt/meter)");
title("Electric Field Plot (Abrupt PN Junction)");

figure(3);
plot(dist,[Doping(1:round(x_num/2))+p(1:round(x_num/2)); Doping(round(x_num/2)+1:x_num)-n(round(x_num/2)+1:x_num)].*1e-6,"LineWidth",2);
xlabel("Distance (along X axis)");
ylabel("Charge Density (cm-3)");
title("Charge Density Plot (Abrupt PN Junction)");

figure(4);
plot(dist,p.*1e-6,dist,n.*1e-6,"LineWidth",2);  %1e-6 is used to convert m-3 to cm-3
xlabel("Distance (along X axis)");
ylabel("Electron & Hole Concentration (cm-3)");
title("Electron, Hole Density Plot (Abrupt PN Junction)");
legend('Hole','Electron',Location='east');



% Energy Band Calculation

% Fermi Level of Semiconductor (Ef) assume to be at 0eV

Ef = zeros(x_num,1); % Fermi Level of Device

Efi_p = Vth * log((p)/ni); % Intrinsic Fermi Level on P Side
Ecb1 = Efi_p + E_g/2;  % Conduction Band Level on P Side
Evb1 = Efi_p - E_g/2;  % Valance Band Level on P Side

Efin = -Vth * log((n)/ni);  % Intrinsic Fermi Level on N Side
Ecb2 = Efin + E_g/2;  % Conduction Band Level on N Side
Evb2 = Efin - E_g/2;  % Valance Band Level on N Side

figure(5);
plot(dist,Ef,"LineWidth",2);
hold on;
plot(dist,[Efi_p(1:round(x_num/2)); Efin(round(x_num/2)+1:x_num)],"LineWidth",2);
plot(dist,[Ecb1(1:round(x_num/2)); Ecb2(round(x_num/2)+1:x_num)],"LineWidth",2);
plot(dist,[Evb1(1:round(x_num/2)); Evb2(round(x_num/2)+1:x_num)],"LineWidth",2);
legend('E Fermi','E Fi','Ec','Ev');
xlabel("Distance (along X axis)");
ylabel("Energy (eV)");
title("Energy Band Diagram (Abrupt PN Junction)");

% Comparison Between Abrupt Junction with Depletion Approximation and Without Depletion Approximation
% Comparision with Assignment 2 (Abrupt Junction)
% Below Code Taken from Assignment - 2
Vbi = Vth * log(Na * Nd / (ni^2));

% Depletion Region Width Calculation
W_dep = sqrt(2 * epsilonsi * (1/Na + 1/Nd) * Vbi / q);
W_n = (W_dep * Na) / (Na + Nd); %Depletion width on n side using depletion approximation
W_p = (W_dep * Nd) / (Na + Nd); %Depletion width on p side using depletion approximation

figure(6);
plot(dist,[Doping(1:round(x_num/2))+p(1:round(x_num/2)); Doping(round(x_num/2)+1:x_num)-n(round(x_num/2)+1:x_num)].*1e-6);
hold on;

% Charge profile
Doping = zeros(1, x_num);
Doping(round(x_num/2)-round(W_p/dx)+1:round(x_num/2)) = -Na;
Doping(round(x_num/2)+1:round(x_num/2)+1+round(W_n/dx)) = Nd;

plot(dist,Doping.*1e-6,"LineWidth",2);
legend('Without Depletion Approximation','With Depletion Approximation',Location='northwest');
xlabel("Distance (along X axis)");
ylabel("Charge Density (cm-3)");
title("Charge Profile Comparison (Abrupt Junction)");

% Poission's Equation RHS
pv = transpose(-(q .* Doping) .* (dx^2 / epsilonsi));

% System of Linear Equations AV = B
A = zeros(x_num, x_num);
A(1,1) = 1;
A(x_num,x_num) = 1;
for i = 2:x_num-1
    A(i,i) = -2;
    A(i,i-1) = 1;
    A(i,i+1) = 1;
end

B = pv;
B(1,1) = -Vbi/2;
B(x_num,1) = Vbi/2;

Vbi_simulated = A\B;

figure(7);
plot(dist,V,"LineWidth",2);
hold on;
plot(dist,Vbi_simulated,"LineWidth",2);
legend('Without Depletion Approximation','With Depletion Approximation',Location='northwest');
xlabel("Distance (along X axis)");
ylabel("Potential (Volt)");
title("Potential Profile Comparison (Abrupt Junction)");
