*RC-Check
********Include Files*******
.OPTION POST=2
.hdl “Resistor.va"
.hdl “Capacitor.va"
.model cap simple_capacitor
.model res simple_resistor
******* Netlist ************
X1 1 2 res resistance=2000 
X2 2 0 cap C=0.0796e-6
Vs 1 0 AC 1V
*************************
.AC DEC 10 1 1MEG
.print V(2) V(1)
.plot V(2) vs V(1)
.end
