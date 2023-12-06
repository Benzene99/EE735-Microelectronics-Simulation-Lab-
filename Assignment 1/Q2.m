clc; clearvars; %clearing all previous variables

X = 1000;Y = 100; %declaring number of pts on the X and Y axes 
m = Y;
n = X;
VPlate1 = -5; %assigning plate voltage
VPlate2 = 5;
Lp1 = 500; %length of plate1
Dp = 6;



for VPlate3 = -5:5:5 %changing plate3 voltage from -5V to +5V
    for D23 = 30:30:90 %changing distance between plates p2 and p3
        Vt = zeros(m,n) ; %initialising matrices for voltage and electric field calculation
        E1= zeros(m,n) ;

        %assigning plate voltages across different plates

        Vt(m/2+(Dp/2-1):m/2+(Dp/2+1),n/2-Lp1/2:n/2+Lp1/2) = VPlate1 ;
        Vt(m/2-(Dp/2-1):m/2-(Dp/2+1),n/2-Lp1/2:n/2-(D23/2)) = VPlate2 ;
        Vt(m/2-(Dp/2-1):m/2-(Dp/2+1),n/2+(D23/2):n/2+Lp1/2)= VPlate3;
        
        
        
        for iter = 1:1:1000
        for i=3:m-2
            for j=3:n-2
                % keeping the plate voltages intact
                %condition for plate3 has to be modelled differently due to
                %the case where Vplate3=0V
                
                if Vt(i,j)==VPlate1 || Vt(i,j)==VPlate2 || (i>=m/2-(Dp/2-1))&&((i>=m/2-(Dp/2+1))&&(j>=n/2+(D23/2))&&(j<=n/2+Lp1/2))
                    break;
                else
                 % Laplace equation assuming equal spacing of points in both axes
                Vt(i,j)=(Vt(i+1,j)+Vt(i-1,j)+Vt(i,j+1)+Vt(i,j-1))/4;
                end
            end
        end
        end
        

        [Ex,Ey]= gradient(Vt); %Finding out electric field using gradient function
        Ex=-Ex;
        Ey=-Ey;
        
     

        E = sqrt(Ex.^2+Ey.^2); %calculating magnitude of electric field
        M= max(E);
        epslion = 8.85*10^-12;
        %initialising charge matrices across the plate dimensions
        delQ2 = zeros(1,length(n/2-Lp1/2:n/2-(D23/2))) ;
        delQ3 = zeros(1,length(n/2+(D23/2):n/2+Lp1/2)) ;
        %calculating charge per unit length using Gauss Law
        for i = m/2-(Dp/2-1)
        for j = n/2-Lp1/2:n/2-(D23/2)
            delQ2(1,j) = E(i,j)*epslion;
        end
        for k = n/2+(D23/2):n/2+Lp1/2
            delQ3(1,k) = E(i,k)*epslion;
        end
        end

        Q22=sum(delQ2);
        Q33=sum(delQ3);
        %calculating capacitance per unit lenth
        C2 = abs(Q22/(VPlate2-VPlate1));
        C3 = abs(Q33/(VPlate3-VPlate1));

        %calculating capacitances for different voltages which result in
        %different configurations
        if VPlate3==-5 
            Ceq=C2;
        elseif VPlate3==0
            Ceq=(C2*C3)/(C2+C3);
        elseif VPlate3==5 
            Ceq=(C2+C3);
        end

        fprintf("The Numerical Value of Capacitance is %d \n",Ceq);
        
        for i=2:m-1
        for j=2:n-1
            if M== sqrt(Ex.^2+Ey.^2)
               fprintf('The maximum value of Electric field occurs at %d and %d \n',i,j);

            end
        end
        end

        figure
        contour(Vt);
        quiver(Ey,Ex,4), hold off ;
        title(' Electrical field of the parallel plate capacitor');

        figure
        pcolor(Vt);
        shading interp; colorbar;
        xlabel(' Lx');
        ylabel(' Ly ');
        title (' Electric potential of the parralel plate capacitor');

        figure
        contour(Vt);
        title(' Equipotential surfaces of the parallel plate capacitor');
    end
end




    


