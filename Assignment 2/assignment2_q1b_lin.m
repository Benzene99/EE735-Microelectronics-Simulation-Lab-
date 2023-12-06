clc , clear all ;
close all;

L = 0.004; 
X = 800;
Y = 200;
ni = 10^10;
Vbi = 0.6 ;
Vt = 0.026 ;
Nd = 1*10^15 ;
Na =  ((ni^2)/Nd)*exp(Vbi/Vt);
q = 1.6*10^-19;
epslon = 11.8*8.85*10^-14 ;

Nedge = Nd;
Wdep = sqrt((6*epslon*Vbi)/(q*Nedge));
Wp = X/2-floor((Wdep/2)*10^6);
Wn = X/2+floor((Wdep/2)*10^6);

DepWidth= Wp-1:Wn+1;
ChargeDensity = zeros(1,X) ;
ElectricField = zeros(1,X);
ElectricPotential = zeros(1,X);

for i=Wp:Wn
        ChargeDensity(1,i)=(q*2*Nedge*((i-X/2)))/(Wn-Wp);
end


DiffMat = zeros(X);
for i=1:X
    for j=1:X
        if i==1&&j==1
            DiffMat(i,j)=1;
        elseif i==X&&j==X
            DiffMat(i,j)=1;
        elseif i==j
            DiffMat(i,j)=-2;
        elseif (i==1&&j==2)||(i==X&&j==X-1)
            DiffMat(i,j)=0 ;
        elseif i==j-1
            DiffMat(i,j)=1;
        elseif i==j+1
            DiffMat(i,j)=1;
        end
    end
end
a=DiffMat;
b=((-q*1/epslon)*ChargeDensity); 

% Diagonalization
for k = 1 : X-1
for i = k+1 : X
a(i,k) = a(i,k)/a(k,k);
for j = k+1:X
a(i,j) = a(i,j)-a(i,k)*a(k,j);
end
end
end
% Forward substitution to solve Ld=b
x(1)=b(1) ;
for i = 2 : X
s=0;
for j = 1 : i-1
s = s + a(i,j)*x(j);
end
x(i) = b(i)-s;
end
% back substitution to solve Ux=d
x(X) = x(X) / a(X, X) ;
for i = X-1 : -1 : 1
s=0 ;
for j = i+1:X
s = s + a( i, j) * x(j) ;
end
x(i) = (x(i) - s)/a(i,i) ;
end

ElectricPotential = x';
subplot(2,1,1);
plot(1:X,ChargeDensity);
xlabel('Distance');
ylabel('Charge Density') ;
title = ('ChargeDensity Vs Distance') ;
grid;

subplot(2,1,2);
plot(1:X,ElectricPotential);
xlabel('Distance');
ylabel('ElectricPotential') ;
title = ('ElectricPotential Vs Distance') ;
grid ;