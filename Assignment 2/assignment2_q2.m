clc ;
clear vars ;
clear all;
%Given parameters
Nd= 1e15;
Na= 1e17;

%Assumptions made to calculate Fermi Level
Vt = 0.0258; % Thermal voltage at room temperature
Ev = 0;% uppermost level of valence band assumed to be at 0 eV
Ea= Ev + 0.067; %Aluminium acceptor ionization energy level in silicon 
Ec = 1.12; %silicon band gap is 1.12 eV
Ed = Ec - 0.054; % Arsenic donor ionization energy level in silicon 

% Source for ionization energy levels- https://inst.eecs.berkeley.edu/~ee130/sp07/lectures/lecture3.ppt

Nc = 3e19; %Effective Density of States in the Conduction Band (NC) per cm^-3
Nv = 1e19; %Effective Density of States in the Valence Band (NV) per cm^-3

Efi = (Ec+Ev)/2 ; % intrinsic fermi level

Ef = Ea; 

Ni = 2000;% no of iterations for Newton Raphson method

for z =1:Ni 

 Nd1 = (Nd/(1+2*exp((Ef-Ed)/Vt)));
 Na1 = (Na/(1+4*exp((Ea-Ef)/Vt)));
 n = (Nc*(exp((Ef-Ec)/Vt)));
 p = (Nv*(exp((Ev-Ef)/Vt)));

 q = Nd1 - Na1 + p - n ;

 q1 = -(Nd1*(1+0.5*exp((Ed-Ef)/0.0258))^(-1)*(0.0258)^(-1))-(Na1*(1+0.25*exp((Ef-Ea)/0.0258))^(-1)*(0.0258)^(-1))-(n*(0.0258)^(-1))-(p*(0.0258)^(-1));
 Ef = Ef - (q/q1);

end

fprintf( "The fermi level of the compensated semiconductor is  %f eV above the Valence Band\n",Ef);


X= 1:10;
EC = zeros(1,10);
EV = zeros(1,10);
EFi = zeros(1,10);
EF = zeros(1,10);

 for i=1:10
     EC(1,i) = Ec;
     EV(1,i) = Ev;
     EFi(1,i) = Efi;
     EF(1,i) = Ef;
 end

 plot (X,EC,'Displayname','Ec');
 hold on
 plot (X,EFi,'Displayname','Efi');
 plot (X,EF,'Displayname','Ef');
 E1 = plot (X,EV,"r",'Displayname','Ev');
 h = gca;
 h.XAxis.Visible = 'off';
 title("Energy Band Diagram")
 ylabel("(E - Ev) in eV");
 legend;

 
     



