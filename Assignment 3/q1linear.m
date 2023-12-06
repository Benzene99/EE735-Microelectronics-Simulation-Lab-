% Q1.a- Linear PN Juntion
clc;
clear vars;
clear all;


L = 1.5; % Length of both n and p region in micrometer
epsilonrsi = 11.8; % Dielectric constant of Silicon
epsilono = 8.84e-14; % Permittivity of free space in F/cm
ni = 1.5e10; % Intrinsic carrier concentration per cm^3
Na = 4e15; % Acceptor doping concentration per cm^3
Nd = 4e15; % Donor doping concentration per cm^3
T = 300; % Temperature in K
q = 1.6e-19; % Charge of an electron
Vth = 0.026; % Thermal voltage in volts
lingrad = ((Nd+Na)/L)*10^4; % Impurity gradient (per cm^4) in linerly doped region
Vbi_calculated = Vth*log((Na*Nd)/ni^2); % Using Boltzmann Law
W = ((12*epsilonrsi*epsilono*Vbi_calculated)/(q*lingrad))^(1/3); % Depletion width in cm
W_n = 10^4*0.5*W; % Depletion width along n-side in micrometer taken till W/2 as mentioned
W_p = 10^4*0.5*W; % Depletion width along p-side in micrometer taken till W/2 as mentioned

% Depletion Approximation Case

dx = 1e-3; % Grid step 
X_num = -L:dx:L; % X-axis

% Volume Charge Density matrix in C/cm^3
rho = zeros(1,length(X_num));
rho(1,round((-W_p+L)/dx+1):round((W_n+L)/dx+1)) = q*lingrad*X_num(1,round((-W_p+L)/dx+1):round((W_n+L)/dx+1))*10^(-4);

% Electric Field Profile
E_calc = zeros(1,length(X_num));

for i = 2:length(X_num)
    E_calc(1,i) = E_calc(1,i-1) + dx * 0.5 * (rho(1,i-1) + rho(1,i));
end

E_calc = 10^(-4) * (epsilonrsi * epsilono)^(-1) * E_calc;

% Voltage Profile
V_calc = zeros(1,length(X_num));

for i = 2:length(X_num)
    V_calc(1,i) = V_calc(1,i-1) + dx * 0.5 * (E_calc(1,i-1) + E_calc(1,i));
end

V_calc = -(10^(-4))*V_calc;

% Without Depletion Approximation Case

h = 1e-2; % Step size 
x_num2 = -L:h:L; % X-axis

% Potential matrix
V = zeros(length(x_num2),1);

% Charge Density matrix assuming complete ionization
rho_v = zeros(length(x_num2),1);

% Generating A matrix
A = zeros(length(x_num2),length(x_num2));
A(1,1) = 1;
A(1,2) = -1;
A(length(x_num2),length(x_num2)) = 1;
A(length(x_num2),length(x_num2)-1) = -1;
for i = 2:length(x_num2)-1
    A(i,i-1) = 1;
    A(i,i) = -2;
    A(i,i+1) = 1;
end

% Matrix for dB/dV
B1 = zeros(length(x_num2),length(x_num2));

for z = 1:400 

    V(1,1) = -Vth*log(Na/ni); % Boundary Condition for p-region end
    V(length(x_num2),1) = Vth*log(Nd/ni); % Boundary Condition for n-region end
    
    p = ni * exp(-V/Vth);
    n = ni * exp(V/Vth);

    % Charge Density Matrix
    rho_v(1,1) = 0; % Boundary Condition for p-region end
    rho_v(2:round((0.5*L)/h),1) = q*(p(2:round((0.5*L)/h),1) - n(2:round((0.5*L)/h),1) - Na);
    rho_v(round((0.5*L)/h+1):round((1.5*L)/h+1),1) = q*(p(round((0.5*L)/h+1):round((1.5*L)/h+1),1) - n(round((0.5*L)/h+1):round((1.5*L)/h+1),1) + lingrad*x_num2(1,round((0.5*L)/h+1):round((1.5*L)/h+1))'*10^(-4));
    rho_v(round((1.5*L)/h+2):length(x_num2)-1,1) = q*(p(round((1.5*L)/h+2):length(x_num2)-1,1) - n(round((1.5*L)/h+2):length(x_num2)-1,1) + Nd);
    rho_v(length(x_num2),1) = 0; % Boundary Condition for n-region end
    
    % Newton Raphson using Jacobian matrix
    B = -(h^2*10^(-8)*(epsilonrsi*epsilono)^(-1))*rho_v;
    f = A * V - B;
    
    for i = 1:length(x_num2)
        B1(i,i) = -q*(h^2*10^(-8)*(epsilonrsi*epsilono)^(-1))*(Vth)^(-1)*ni*(-exp(-V(i,1)/Vth)-exp(V(i,1)/Vth));
    end

    J = A - B1;

    % LU Decomposition to calculate inverse( to bypass matrix singularity
    % error)
    LT = eye(length(x_num2));
    UT = J;
    for i = 2:length(x_num2)
        for j = i:length(x_num2)
            LT(j,i-1) = UT(j,i-1)/UT(i-1,i-1);
            UT(j,:) = UT(j,:) - (UT(j,i-1)/UT(i-1,i-1))*UT(i-1,:);
        end
    end

    y = zeros(length(x_num2),1);
    for i = 1:length(x_num2)
        temp = 0;
        for j = 1:length(x_num2)
            if(i~=j)
                temp = temp + LT(i,j)*y(j,1);
            end
        end
        y(i,1) = f(i,1) - temp;
    end
    
    x = zeros(length(x_num2),1);
    for i = length(x_num2):-1:1
        temp = 0;
        for j = 1:length(x_num2)
            if(i~=j)
                temp = temp + UT(i,j)*x(j,1);
            end
        end
        x(i,1) = (y(i,1) - temp)/UT(i,i);
    end

    V = V - x;

end

% Charge Density 
p = ni * exp(-V/Vth);
n = ni * exp(V/Vth);
rho_v(1:round((0.5*L)/h),1) = q*(p(1:round((0.5*L)/h),1) - n(1:round((0.5*L)/h),1) - Na);
rho_v(round((0.5*L)/h+1):round((1.5*L)/h+1),1) = q*(p(round((0.5*L)/h+1):round((1.5*L)/h+1),1) - n(round((0.5*L)/h+1):round((1.5*L)/h+1),1) + lingrad*x_num2(1,round((0.5*L)/h+1):round((1.5*L)/h+1),1)'*10^(-4));
rho_v(round((1.5*L)/h+2):length(x_num2),1) = q*(p(round((1.5*L)/h+2):length(x_num2),1) - n(round((1.5*L)/h+2):length(x_num2),1) + Nd);
rho_v = rho_v';

% Potential in Volts
V = V';

% Electric Field Intensity in Volts/cm
E = zeros(1,length(x_num2));
for i = 2:length(x_num2)
    E(1,i) = E(1,i-1) + h * 0.5 * (rho_v(1,i-1) + rho_v(1,i));
end
E = 10^(-4) * (epsilonrsi * epsilono)^(-1) * E;

% Potential, Electric Field, Charge Concentration Profile

figure(1);
plot(X_num* 10^(-4),V_calc,'Displayname','with approx.',"LineWidth",2);
hold on;
plot(x_num2* 10^(-4),V-V(1,1),'Displayname','without approx.',"LineWidth",2);
xlabel('Distance ( along X axis )');
ylabel('V (Volts)');
title('Voltage Profile');
legend;
grid;
hold off;

figure(2);
plot(X_num* 10^(-4),E_calc,'Displayname','with approx.',"LineWidth",2);
hold on;
plot(x_num2* 10^(-4),E,'Displayname','without approx.',"LineWidth",2);
xlabel('Distance ( along X axis )');
ylabel('E (Volts/cm)');
title('Electric Field Profile');
legend;
grid;
hold off;

figure(3);
plot(X_num* 10^(-4),rho/q,'Displayname','with approx.',"LineWidth",2);
hold on
plot(x_num2* 10^(-4),rho_v/q,'Displayname','without approx.',"LineWidth",2);
xlabel('x ( cm )');
ylabel('Charge Density(cm-3)');
title('Charge Density Profile');
legend;
grid;
hold off;


% Electron and Hole concentration

p = ni * exp(-V/Vth);
n = ni * exp(V/Vth);

figure(4);
semilogy(x_num2* 10^(-4),p,'Displayname','hole',"LineWidth",2);
hold on;
semilogy(x_num2* 10^(-4),n,'Displayname','electron',"LineWidth",2);
xlabel('x (along X axis )');
ylabel('Carrier Density (cm-3 )');
title('Electron and Hole Densities');
legend;
grid;
hold off;

% Energy Band Diagram Plot

E_F = zeros(1,length(x_num2));
E_i = E_F - V;
E_g = 1.12; 
E_C = E_i + E_g/2;
E_V = E_i - E_g/2;

figure(5);
plot(x_num2,E_C,'Displayname','E_C',"LineWidth",2);
hold on
plot(x_num2,E_i,"--",'Displayname','E_{mid}',"LineWidth",2);
plot(x_num2,E_F,"k",'Displayname','E_F',"LineWidth",2);
plot(x_num2,E_V,"r",'Displayname','E_V',"LineWidth",2);
h = gca;
h.XAxis.Visible = 'off';
ylabel("(E-Ef) in eV");
title("Energy Band Diagram");
legend;
grid;
hold off;