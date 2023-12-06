%Assignment 4
%Q1a

clc;
clear vars;

L= 10e-4; %defined a region of length 10um in cm
dx=0.01e-4; %step size for calculations
dist=0:dx:L; %distance vector along X-axis
gridno= round(L/dx)+1; %no of steps in our grid depending on step size
A = zeros(gridno,gridno); %initialising matrix to calculate discrete double derivative 
A(1,1) = 1; %to give us proper boundary condition of N matrix at first node
A(gridno,gridno) = 1; %to give us proper boundary condition of N matrix at last node
tau=1e-7; %as given in question
D=0.1; %as given in question


%Generating the A matrix for numerical solution on differential equation
for i= 2:gridno-1
    
    A(i,i-1) = (1/dx^2);
    A(i,i)= -((2/dx^2) +1/(D*tau));
    A(i,i+1)= (1/dx^2);
end

%Initialsing N,B1,B2 matrices for particle and two given boundary
%conditions
N1 = zeros(gridno,1);
N2= zeros(gridno,1);
B1 =zeros(gridno,1);
B2 =zeros(gridno,1);
B1(1,1) = 1e12;%to give us proper boundary condition of N matrix at first node
B1(gridno,1)=0;%first boundary condition
B2(1,1) = 1e12;%to give us proper boundary condition of N matrix at first node
B2(gridno,1)=1e12*exp(-10);%boundary condition at B for second condition by solving differential equation of J=-Ddn/dx and putting value of 10um

N1 = A\B1;
J1 = -D*(gradient(N1))/dx;

%plots for 2 boundary conditions in absolute and logarithmic scale
figure(1);
plot(dist,N1,"LineWidth",2);
xlabel("Distance (in cm)");
ylabel("Number of particles (per cm-3)");
title("Diffusion of particles for first Boundary Condition");

figure(2);
plot(dist,J1,"LineWidth",2);
xlabel("Distance (in cm)");
ylabel("Flux of particles (per cm-3 s-1)");
title("Flux of particles for first boundary condition");


N2 = A\B2;
J2 = -D*(gradient(N2))/dx;
p=J2(gridno);
fprintf('The particle flux at B for second boundary condition is %d cm^-^2s^-^1',p);



figure(3);
plot(dist,N2,"LineWidth",2);
xlabel("Distance (in cm)");
ylabel("Number of particles (per cm-3)");
title("Diffusion of particles for second Boundary Condition");

figure(4);
plot(dist,J2,"LineWidth",2);
xlabel("Distance (in cm)");
ylabel("Flux of particles (per cm-3 s-1)");
title("Flux of particles for second boundary condition");

figure(5)
plot(dist,log10(N1),'LineWidth',2,'LineStyle','-');
hold on
plot(dist,log10(N2),'g','LineWidth',2,'LineStyle','--');
hold off
grid on;
title('Particle Profile Semilog scale');
xlabel('x (cm) ');
ylabel('log(Particle Concentration (cm^-^3)) ');
legend('Boundary Condition at B: n=0','Boundary Condition at B: J=kn');

figure(6)
plot(dist,log10(J1),'LineWidth',2,'LineStyle','-');
hold on
plot(dist,log10(J2),'g','LineWidth',2,'LineStyle','--');
hold off
grid on;
title('Particle Flux Semilog scale');
xlabel('x (cm) ');
ylabel('log(Particle Flux (cm^-^2s^-^1)) ');
legend('Boundary Condition at B: n=0','Boundary Condition at B: J=kn');
