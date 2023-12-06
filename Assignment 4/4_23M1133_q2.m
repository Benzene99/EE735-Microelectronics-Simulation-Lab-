%Assignment 4
%Q2

clc;
clear vars;
L = 100e-4; % defined a region of length 100um in cm
D = 1e-4; % Diffusion constant as given
N0 = 2000; % Particle density at x=0,t=0 as given
dx = 1e-4; %%step size for calculations
C = 0.7; % Constant for Euler method
ti = 1e-3;% time gap of 1ms taken from impulse injection
tf = 15e-3;%final time is taken as 15 ms for accurate results
dt = (C*dx^2)/D;%as derived in PPT
dist = 0:dx:L; %distance vector
T = 0:dt:tf;%time vector
Ni=5000;%no of iterations
T1=length(T);
dist1=length(dist);

%initialising particle matrix for each time instant as a column
N = zeros(dist1,T1);
N(1,:) = N0;%Initialising first row at every time instant as N0

%Loop iterations to find N(i,t) as derived in PPT
for z = 1:Ni
    for t = 2:T1
        for i = 2:dist1-1
            N(i,t) = (C*(N(i-1,t)+N(i+1,t))+N(i,t-1))/(2*C+1);
        end
    end
end

%Plots of particle evolution in space and time
figure;
plot(dist,N(:,floor(ti/dt+1)),'Displayname','t = 1');
hold on;
plot(dist,N(:,floor(2*ti/dt+1)),'Displayname','t = 2');
plot(dist,N(:,floor(3*ti/dt+1)),'Displayname','t = 3');
plot(dist,N(:,floor(4*ti/dt+1)),'Displayname','t = 4');
plot(dist,N(:,floor(5*ti/dt+1)),'Displayname','t = 5');
plot(dist,N(:,floor(6*ti/dt+1)),'Displayname','t = 6');
plot(dist,N(:,floor(7*ti/dt+1)),'Displayname','t = 7');
plot(dist,N(:,floor(8*ti/dt+1)),'Displayname','t = 8');
plot(dist,N(:,floor(9*ti/dt+1)),'Displayname','t = 9');
plot(dist,N(:,floor(10*ti/dt+1)),'Displayname','t = 10');
plot(dist,N(:,floor(11*ti/dt+1)),'Displayname','t = 11');
plot(dist,N(:,floor(12*ti/dt+1)),'Displayname','t = 12');
plot(dist,N(:,floor(13*ti/dt+1)),'Displayname','t = 13');
plot(dist,N(:,floor(14*ti/dt+1)),'Displayname','t = 14');
plot(dist,N(:,floor(15*ti/dt+1)),'Displayname','t = 15');
title("Particle Profile Using Numerical Solution");
xlabel("x in cm ");
ylabel("Concentration per cm^-3 ");
legend;
grid on;
hold off;

 
% Analytical solution
Analytical_n = (N0) * erfc(dist / (2 * sqrt(D * ti)));

figure;
plot(dist, Analytical_n, 'DisplayName', 'Analytical Solution');
hold on;
plot(dist, N(:, floor(ti/dt+1)), 'DisplayName', 'Numerical Solution');
xlabel('Length (cm)');
ylabel('Particle Concentration (per cm^3)');
title('Comparison between Numerical and Analytical solution');
legend;
grid on;
hold off;