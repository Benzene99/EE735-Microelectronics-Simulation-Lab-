clc; 
clear vars ;
close all;
%declaring boundary of device
X = 800;
Y = 200;

%declaring given parameters
Na = 7.48e+14 ;
Nd = 1e+15 ;
Vbi = 0.6 ;
q = 1.6e-19;
epsilonSi = 11.8*8.85*10^(-14) ;

%calculation of depletion widths
Wd = sqrt((4*epsilonSi*Vbi)/(q*Na)); 
Wp = X/2-floor((Wd/2)*10^5);
Wn = X/2+floor((Wd/2)*10^5);
 
DepletionWidth= Wp:Wn;
m = X ;
rho = zeros(1,m) ; %%initialising charge,potential,field matrix
Field = zeros(1,m);
Potential = zeros(1,m);

rho(Wp:X/2)=-Na;
rho(X/2+1:Wn)=Nd;

%Initialising the finite difference Matrix

FinDiff = zeros(m);
for i=1:m
    for j=1:m
        if i==1&&j==1
            FinDiff(i,j)=1;
        elseif i==m&&j==m
            FinDiff(i,j)=1;
        elseif i==j
            FinDiff(i,j)=-2;
        elseif (i==1&&j==2)&&(i==m&&j==m-1)
            FinDiff(i,j)=0 ;
        elseif i==j-1
            FinDiff(i,j)=1;
        elseif i==j+1
            FinDiff(i,j)=1;
        end
    end
end
a=FinDiff;
b=(((q*10^-10)/epsilonSi)*rho)';

%Performing Gauss Siedel Method

x=Potential';
%initialising Lower triangular,Upper triangular and diagonal matrices
L = zeros(m);
U = zeros(m);
D = zeros(m);

for i=1:m
    for j=1:m
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