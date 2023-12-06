%Assignment 4
%Q1b

clc;
clear vars;

L= 10e-4;%defined a region of length 10um in cm
dx=0.01e-4;%step size for calculations
dist=0:dx:L;%distance vector along X-axis
gridno= round(L/dx)+1;%no of steps in our grid depending on step size
A = zeros(gridno,gridno);%initialising matrix to calculate discrete double derivative 
A(1,1) = 1;%to give us proper boundary condition of N matrix at first node
A(gridno,gridno) = 1;%to give us proper boundary condition of N matrix at last node
tau=1e-7;%as given in question
D=0.1;%as given in question
J1=1e13; %Particle flux as given
xflux = 3 + 0.5*3; % Point of Particle Flux;Last digit of roll no = 3
xflux = xflux*10^-4;  % convert to cm
fluxpt = round(1 + xflux/dx); %No of points upto the point of flux introduction

%Generating the A matrix for numerical solution on differential equation
for i= 2:gridno-1
    
    A(i,i-1) = (1/dx^2);
    A(i,i)= -((2/dx^2) +1/(D*tau));
    A(i,i+1)= (1/dx^2);
end

%to give us proper boundary condition at the point of flux introduction
A(fluxpt,fluxpt) = 1;
A(fluxpt,fluxpt-1) = -1;
A(fluxpt,fluxpt+1) = 0;

%Initialising N,B matrices for particles
N = zeros(gridno,1);
B =zeros(gridno,1);
%boundary conditions as given in question
B(1,1) = 0;
B(gridno,1) = 0;
B(fluxpt,1) = (dx*J1)/D;

%solving for N by inverse of A matrix
N = A\B;
J = -D*(gradient(N))/dx;
p1= J(1);
p2= J(gridno);

fprintf('The particle flux at A is %d cm^-^2s^-^1\n',p1);
fprintf('The particle flux at B is %d cm^-^2s^-^1',p2);


%Plots for particle profile and flux
figure(1);
plot(dist,N,"LineWidth",2);
grid on;
xlabel("Distance (in cm)");
ylabel("Number of particles (per cm-3)");
title("Diffusion of particles with given Boundary Conditions");

figure(2);
plot(dist,J,"LineWidth",2);
grid on;
xlabel("Distance (in cm)");
ylabel("Flux of particles (per cm-3 s-1)");
title("Flux of particles with given boundary conditions");

%Analytical Solution

N1 = 1e10;

%Left side of gd_pos 0 < x < 4e-4
N2 = -(N1/(exp(-4)+exp(4)));
N3 = N1/(exp(-4)+exp(4));
An_n(1,1:gd_pos) = N2*exp(-x(1,1:gd_pos)/l)+N3*exp(x(1,1:gd_pos)/l);

%Right side of gd_pos 4e-4 < x < 10e-4
C3 = N1/(1+exp(-12));
C4 = -(N1/(1+exp(12)));
An_n(1,gd_pos:length(x)) = C3*exp(-(x(1,gd_pos:length(x))-pos_x)/l)+C4*exp((x(1,gd_pos:length(x))-pos_x)/l);

figure(3)
plot(x, n, 'b', x, An_n, 'r',"LineWidth",2);
xlabel('Length(in cm)---------->');
ylabel('Particle Concentration(n)(in/cm^3)---------->')
title('Comparison of Numerical and Analytical Solutions');
legend('Numerical Solution','Analytical Solution');