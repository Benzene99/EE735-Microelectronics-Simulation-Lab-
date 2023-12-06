* Basic mos_test
******** Include Files
.OPTION POST=2 INGOLD=2
.hdl "NMOS.va" $ change the names according to your file_names
.model mos simple_nmos $ change module name accordingly 
****NETLIST
X1 3 2 0 0 mos
Vdum 4 3 DC 0 
Vdd 4 0 DC 2 
*Reference Hspice code for NMOS
Vgg 2 0 DC 2
*****Analysis $ LV HV td tr tf PW PER
.dc Vdd 0 2 0.02 vgg 0.6 2 0.35 $ to plot Output Characteristics
*.dc Vgg 0 2 0.02 vdd 0 2 0.5 $ to plot Input Characteristics
*For Ip char uncomment line 14 and for Op char uncomment line 13 
.print I(Vdum) 
