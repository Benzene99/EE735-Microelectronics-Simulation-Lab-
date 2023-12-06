%Assignment 4
%Q3

q = 1.6e-19;              % Charge of an electron (Coulombs)
N_D = 1e15;               % Donor concentration (per cm^3)
N_A = 1e15;               % Acceptor concentration (per cm^3)
L = 0.3;                  % Length of the diode (cm)
Wd = L / 2;
up = 450;                 % Hole mobility (cm^2/Vs)
taup = 20e-6;             % Hole lifetime (s)
Ni = 1.5e10;              % Intrinsic carrier concentration (per cm^3)
Jr = 2000;                % Current density (A/cm^2)
gridsize = 100;           % Number of grid points in the N-side of the diode
instants = 30000;         % Number of time instants
Po = (Ni^2) / N_A;        % Hole concentration at equilibrium without applied voltage
er = 11.8;
eo = 8.854e-14;
esi = er * eo;
Va = 1;                    % Applied voltage (Volts)
Vt = 0.026;                %Thermal Voltage
dist = 0;                  
gridpoints = linspace(0, 1e-3, gridsize);

% Calculating depletion width in the N-side
Vbi = Vt * log(N_D * N_A / (Ni^2));
Wdep = sqrt((2 * esi / q) * (1 / N_D + 1 / N_A) * (Va + Vbi));
Wn = N_D * Wdep / (N_D + N_A);

Dp = up * Vt; %Using einstein relationship
Ldiff = sqrt(Dp * taup); %Diffusion length calculation
gridstep = (Wd - Wn) / gridsize;

xpositions = linspace(0, Wd, gridsize); 

% Calculating storage time analytically
ts = taup * log((Jr + Jr) / Jr);

% Calculating time step
delt = 15 * ts / instants;

% Initializing matrices
p1 = zeros(gridsize, instants);
delp1 = zeros(gridsize, instants);
delQ = zeros(1, instants);

% Initializing hole concentration matrix at equilibrium
p1(:, :) = Po;

% Applying boundary condition at the depletion edge
delp1(2, 1) = gridstep * Jr / (q * Dp) + delp1(1, 1);

% Initializing delta voltage
for i = 1:gridsize
    delp1(i, 1) = Po * exp(1/Vt) * exp(-dist / Ldiff);
    dist = dist + gridstep;
end

% Running time iterations
K1 = 1 / delt + 2 * Dp / gridstep^2 + 1 / taup;
K2 = gridstep * Jr / (q * Dp);

for j = 2:instants
    for iter = 1:100
        for i = 2:gridsize-1
            delp1(i, j) = (Dp * delp1(i+1, j) / gridstep^2 + ...
                                              Dp * delp1(i-1, j) / gridstep^2 + ...
                                              delp1(i, j-1) / delt) / K1;
        end
    end
end

% Applying boundary condition at the depletion edge
for j = 2:instants
    delp1(1, j) = delp1(2, j) - K2;
end

% Plotting
figure(10);
plot(gridpoints, delp1(:, 1));
title('Hole Concentration Evolution in N-side of Diode');
xlabel('Distance (cm)');
ylabel('Concentration (per cm^3)');
hold on;

% Plotting at selected time instants up to the calculated storage time
for fig = 2:100:instants
    plot(gridpoints, delp1(:, fig));
end

hold off;
