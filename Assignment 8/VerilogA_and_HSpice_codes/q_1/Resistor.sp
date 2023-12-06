* Basic Resistor 
*******************************
* Include Files (* means comment)
*******************************
.OPTION POST=2 INGOLD=2 $ (text after $ is comment)
.hdl “Resistor.va“ $ calling the verilog-a file
.model res simple_resistor $ creating a model instance
*******************************
* Netlist
*******************************
X1 1 0 res resistance=100 $ instantiating the instance as sub-circuit
Vs 1 0 DC 0.8
*******************************
.dc Vs 0 4 0.2
.print V(1)
.print I(Vs)
.end