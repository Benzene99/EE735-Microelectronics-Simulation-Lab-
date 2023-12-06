* Basic mos_test
******** Include Files
.OPTION POST=2 INGOLD=2
.hdl "NMOS.va" $ change the names according to your file_names
.model nmos simple_nmos $ change module name accordingly 
.hdl "PMOS.va" $ change the names according to your file_names
.model pmos simple_pmos $ change module name accordingly 
**************netlist
X1 1 2 3 3 pmos
X2 1 2 0 0 nmos
vdd 3 0 DC 2
c 1 0 0.00050pF
Vgg 2 0 PULSE(0 2 5n 2.5n 2.5n 20n 45n)

.tran 250p 145n
.print V(2) V(1)