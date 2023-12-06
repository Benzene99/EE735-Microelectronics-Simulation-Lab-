clc,clearvars;
xm=1601;
ym=1601;
midx=ceil(xm/2);
midy=ceil(ym/2);
V=zeros(xm,ym); %initialising voltage matrix
Ni=100;   %number of iterations
x=(1:xm)-midx;     %initialising axis about zero
y=(1:ym)-midy;
epsilon=8.85*10^(-12);
Vplate=50;
q=0; %initialising charge variable

%initialising theoretical,simulated and parasitic capacitor arrays
C=zeros(1,31); 
Cth=zeros(1,31);
Cparas=zeros(1,31);
k=1;

plate_1=midy+4; %position of plate below x axis
plate_2=midy-4; %position of plate above x axis


for plate_length=20:50:1520
 for z = 1:Ni    % Loop of iterations
        
        for i=2:xm-1
         for j=2:ym-1     
            
                 V(plate_1,midx-(plate_length/2):midx+(plate_length/2)) = 50; %Fixing the potentials on plates
                 V(plate_2,midx-(plate_length/2):midx+(plate_length/2)) = -50; 
                
                V(i,j)=(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1))/4;
        end
        end
        [Ex,Ey]=gradient(V); %Finding out electric field using gradient function
        Ex=-Ex;
        Ey=-Ey;
        E=sqrt(Ex.^2+Ey.^2); %calculating magnitude of electric field

        
        
end
for j=midx-(plate_length/2):1:midx+(plate_length/2)
             q=q+E(midy+4,j)*epsilon;  %calculating charge on capacitor plate 
 end

        C(k)=q/50;
        Cth(k)=(epsilon*plate_length)/8;
        q=0;
        k=k+1;
end
disp(k);
p=(20:50:1520);
figure(1)
plot(p,C);
title(' Absolute value of capacitor');
hold on
for m=1:31
    Cparas(m)=C(m)-Cth(m);
end
figure(2)
plot(p,Cparas);
title('Value of parasitic capacitance');