* Basic mos_test
******** Include Files
.OPTION POST=2 INGOLD=2
.hdl "NMOS.va" $ change the names according to your file_names
.model mos simple_nmos $ change module name accordingly 
****NETLIST
X1 3 2 0 0 mos
R 4 3 100K
Vdd 4 0 DC 2 

Vgg 2 0 PULSE(0 2.5 5n 2.5n 2.5n 20n 45n)


*****Analysis $ LV HV td tr tf PW PER

.tran 250p 145n

.print V(3) V(2)
