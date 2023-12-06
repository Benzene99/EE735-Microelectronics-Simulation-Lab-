%Assignment 4
%Q1c

clc;
clear vars;

L = 10e-4; % defined a region of length 10um in cm

J = 34e6;%Roll no is 23M1133 so flux is written accordingly
xflux = 5*10^(-4); % Point of Injection
t0 = 10e-6;% time gap of 10us taken from impulse injection
tn = 150 * 10^(-6); %final time is taken as 150 us for accurate results
D = 1e-4; %as given in question
dx = 1e-5; % Step size along X-axis
dist= 0:dx:L; %distance vector
C=0.7; %Implicit condition of C>0.5
dt = (C*dx^2)/D;%as discussed in class
fluxpt = round(1 + xflux/dx); %No of points upto the point of flux introduction
T = 0:dt:tn; %time vector
T1=length(T);
dist1=length(dist);
N = zeros(dist1,T1);%initialising particle matrix for each time instant as a column
N(fluxpt,1) = J;%boundary condition as given
Ni=5000;%no of iterations

%Loop iterations to find N(i,t) as derived in PPT
for z = 1:Ni
    for t = 2:T1
        for i = 2:dist1-1
            N(i,t) = (C*(N(i-1,t)+N(i+1,t))+N(i,t-1))/(2*C+1);
        end
    end
end


%Plot of particle evolution vs time
figure;
plot(dist,N(:,floor(t0/dt+1)),'LineWidth',2,'Displayname','t = 1');
hold on;
plot(dist,N(:,floor(2*t0/dt+1)),'LineWidth',2,'Displayname','t = 2');
plot(dist,N(:,floor(3*t0/dt+1)),'LineWidth',2,'Displayname','t = 3');
plot(dist,N(:,floor(4*t0/dt+1)),'LineWidth',2,'Displayname','t = 4');
plot(dist,N(:,floor(5*t0/dt+1)),'LineWidth',2,'Displayname','t = 5');
plot(dist,N(:,floor(6*t0/dt+1)),'LineWidth',2,'Displayname','t = 6');
plot(dist,N(:,floor(7*t0/dt+1)),'LineWidth',2,'Displayname','t = 7');
plot(dist,N(:,floor(8*t0/dt+1)),'LineWidth',2,'Displayname','t = 8');
plot(dist,N(:,floor(9*t0/dt+1)),'LineWidth',2,'Displayname','t = 9');
plot(dist,N(:,floor(10*t0/dt+1)),'LineWidth',2,'Displayname','t = 10');
plot(dist,N(:,floor(11*t0/dt+1)),'LineWidth',2,'Displayname','t = 11');
plot(dist,N(:,floor(12*t0/dt+1)),'LineWidth',2,'Displayname','t = 12');
plot(dist,N(:,floor(13*t0/dt+1)),'LineWidth',2,'Displayname','t = 13');
plot(dist,N(:,floor(14*t0/dt+1)),'LineWidth',2,'Displayname','t = 14');
plot(dist,N(:,floor(15*t0/dt+1)),'LineWidth',2,'Displayname','t = 15');
title("Particle Profile Using Numerical Solution");
xlabel("x in cm )");
ylabel("Concentration per cm^-3 )");
legend;
grid on;
hold off;

% Analytical solution
analytical_n = J * exp(-(dist - xflux).^2 / (4 * D * dx));
figure(2);
plot(dist, analytical_n, 'b', dist, N(:, floor(t0 / dt + 1)) * 15, 'g');
xlabel('Length (cm)');
ylabel('Particle Concentration (per cm^3)');
title('Comparison between Numerical and Analytical solution');
legend('Analytical Solution', 'Numerical Solution');
grid on;
