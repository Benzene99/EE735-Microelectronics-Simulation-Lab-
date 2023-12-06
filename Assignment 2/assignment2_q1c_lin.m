clc ;
clear vars ;
clear all;

%declaring boundary of device
X = 800;
Y = 200;
%declaring given parameters
ni = 1e+10;
Vbi = 0.6 ;
Vt = 0.026 ;
Nd = 1e+15 ;
Na = 7.48e+14;
q = 1.6e-19;
epslon = 11.8*8.85e-14 ;

%calculation of depletion widths
Nedge = Na;
Wd = sqrt((6*epslon*Vbi)/(q*Nedge));
Wp = X/2-floor((Wd/2)*10^5);
Wn = X/2+floor((Wd/2)*10^5);

DepWidth= Wp-1:Wn+1;
rho = zeros(1,X) ; %initialising charge,potential,field matrix
Field = zeros(1,X);
Potential = zeros(1,X);

for i=Wp:Wn
        rho(1,i)=(q*2*Nedge*((i-X/2)))/(Wn-Wp);
end

%Initialising the finite difference Matrix
FinDiff = zeros(X);
for i=1:X
    for j=1:X
        if i==1&&j==1
            FinDiff(i,j)=1;
        elseif i==X&&j==X
            FinDiff(i,j)=1;
        elseif i==j
            FinDiff(i,j)=-2;
        elseif (i==1&&j==2)&&(i==X&&j==X-1)
            FinDiff(i,j)=0 ;
        elseif i==j-1
            FinDiff(i,j)=1;
        elseif i==j+1
            FinDiff(i,j)=1;
        end
    end
end
a=FinDiff;
x=Potential';
b=((q*10^(9)/epslon)*rho)';

%Performing Gauss Siedel Method
xn=Potential';

%initialising Lower triangular,Upper triangular and diagonal matrices
L = zeros(X);
U = zeros(X);
D = zeros(X);

for i=1:X
    for j=1:X
        if i<j
            L(i,j) = a(i,j);
        elseif i>j
            U(i,j) = a(i,j);
        elseif i==j
            D(i,j) = a(i,j);
        end
    end
end


for i=1:2000 %no of iterations for Gauss Siedel 
    xn = -((inv(D+L))*U*x)+((inv(D+L))*b);
    x=xn;
end

Potential = x';

plot(Wp-1:Wn+1,Potential(Wp-1:Wn+1));
xticklabels([Wp-1:Wn+1]) 
set ( gca, 'xdir', 'reverse' )
xlabel('Lenth along X');
ylabel('ElectricPotential') ;
title = ('ElectricPotential Vs Distance') ;
grid;
