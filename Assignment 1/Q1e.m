clc, clearvars; %clearing all previous variables
X = 1000 ; Y = 60;%declaring number of pts on the X and Y axes 
m = Y;
n = X;
VPlate = 10; %assigning plate voltage
epsilon = 8.85*10^-12;
Ni=1500; %declaring no of iterations
count=1;

Vt = zeros(m,n) ; %initialising matrices for voltage and electric field calculation
Et = zeros(m,n) ;

Vt(m/2-4,n/2-400:n/2+400) = VPlate; %assigning plate voltages across 800nm plate
Vt(m/2+4,n/2-400:n/2+400) = -VPlate ;


for z = 1:Ni
for i=2:m-1
    for j=2:n-1
        if Vt(i,j)==VPlate || Vt(i,j)==-VPlate % keeping the plate voltages intact
            continue;
        elseif (i>=m/2-1)&&(i<=m/2+1)
            Vt(i,j) = (Vt(i-1,j)+5*Vt(i+1,j))/6; %satisfying the necessary boundary conditions midway between the plates
        else
        Vt(i,j)=(Vt(i+1,j)+Vt(i-1,j)+Vt(i,j-1)+Vt(i,j+1))/4;  % Laplace equation assuming equal spacing of points in both axes
        end
    end
end
end



[Ex,Ey]= gradient(Vt); %Finding out electric field using gradient function
Ex=-Ex;
Ey=-Ey;


Et = sqrt(Ex.^2+Ey.^2); %calculating magnitude of electric field
E_plate=Et(m/2-3,n/2-400:n/2+400); %calculating the electric field on the plates


Q = epsilon*sum(E_plate,"all"); %calculating charge per unit length using Gauss Law
C = Q/(2*VPlate); %calculating capacitance per unit lenth
disp(C);

contour(Vt);
hold on, quiver(Ey,Ex,1), hold off ;
title(' Electric field of the parallel plate capacitor');


figure
contour(Vt);
title(' Equipotential surfaces of the parallel plate capacitor');

figure
pcolor(Vt);
shading interp; colorbar;
xlabel(' Lx');
ylabel(' Ly ');
title (' Electric potential of the parralel plate capacitor');

            
    


