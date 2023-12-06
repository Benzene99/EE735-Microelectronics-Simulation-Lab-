clc,clearvars;  %clearing all previous variables
X = 1000;Y = 80; %declaring number of pts on the X and Y axes 
n = X; m = Y;
VPlate=10;%assigning plate voltage
epsilon = 8.85*10^-12;
Ni=1000; %declaring no of iterations
count=1;

V1 = zeros(m,n) ;%initialising matrices for voltage and electric field calculation
E1= zeros(m,n) ;

 V1(m/2-4,n/2-400:n/2+400) = -VPlate ;%assigning plate voltages across 800nm plate
 V1(m/2+4,n/2-400:n/2+400) = VPlate ;

for z = 1:Ni 
        for i=2:m-1
        for j=2:n-1
            
            % keeping the plate voltages intact
            if V1(i,j)==VPlate|| V1(i,j)==-VPlate
                continue;
            end
            % Laplace equation assuming equal spacing of points in both axes
            V1(i,j)=(V1(i+1,j)+V1(i-1,j)+V1(i,j+1)+V1(i,j-1))/4;      
        end
      end
    
    count=count+1;
end
fprintf(' Number of iterations N=%d\n',count);




[Ex,Ey]= gradient(V1); %Finding out electric field using gradient function
Ex=-Ex;
Ey=-Ey;

Et = sqrt(Ex.^2+Ey.^2); %calculating magnitude of electric field

E_plate = Et(m/2-5:m/2-3,n/2-400:n/2+400); %calculating the electric field around the plates cosidering its thickness

Q = epsilon*sum(E_plate,"all"); %calculating charge per unit length using Gauss Law
C = Q/(2*VPlate); %calculating capacitance per unit lenth




disp(C);

contour(V1);
hold on, quiver(Ey,Ex,1), hold off ;
title(' Electric field of the parallel plate capacitor');

figure
contour(V1);
title(' Equipotential surfaces of the parallel plate capacitor');

figure
pcolor(V1);
hold on;
shading interp; colorbar;
xlabel(' Lx');
ylabel(' Ly ');
title (' Electric potential of the parallel plate capacitor');



    


