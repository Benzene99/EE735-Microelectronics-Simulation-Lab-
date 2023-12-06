
Section I: General notes for assignment completion:

1) MOSFET terminals in HSPICE should follow the same order as modeled in Verilog-A(inside the module statement).

2) If using MS Excel for plotting, after copy-paste, use DATA->Text to Columns option to fill the data in columns.

3) To execute circuit(.sp) file, command is hspice -i input.sp -o out Don't give any extension in out. \\
 For errors, see the log in .lis file. Type hspice --help for further details.

4) The print statement with current needs to be put as current from a source. e.g. .print I(VDD)  

5) To output the results in the exponential format use:
.OPTION POST=2 ingold=2

6) In the hspice_sa.pdf manual (link for the same is given in slide 7 of the ppt), \\
the example on JFET compact device model (chapter 31, page 758) might be useful to solve the questions on MOSFET. \\
The output current equation for a MOSFET in linear and saturation regions can be found in any standard textbook on \\
semiconductor device fundamentals, eg., Semiconductor Device Fundamentals by Robert F. Pierret, chapter 17, eqs 17.17 and 17.22. 

7) Try to upload clear and self explanatory figures. Your figures must contain the title and appropriate labels on the X and Y axes.

Please put your graphs and accompanying discussion, observation, etc. in a single .pdf file. The Verilog-A and Hspice codes should be in \\
another file. Zip everything together and submit with the title 'Roll_No_Name_Assignment8.zip'.


Section II: General approach for performing the Hspice simulations:

1) Copy the Verilog-A code from the assignment file given to you and paste in a text editor (eg. notepad++). 
Save it with extension .va. For example save the Verilog-A code for resistor as Resistor1.va. Similarly save the Hspice code for the \\
resistor as Resistor_1.sp.

2) In the Resistor_1.sp file, the Resistor1.va is called through the line: .hdl "Resistor1.va".

3) An instance of the resistor module is created in the Resistor_1.sp file through the line: .model res simple_resistor. \\
Here 'res' is the nickname you give to the instance, and 'simple_resistor' is the module name used in the Resistor1.va file.

4) You can change the value of the parameter 'resistance' either directly in Resistor1.va file or in the Resistor_1.sp file through \\
line: X1 1 0 res resistance=value. Here X1 is a replica of res. You can create many such replicas with different names. The module instance
is treated as sub-circuit in the Hspice file, hence the name of replica always starts with X.

5) Run the Resistor_1.sp file using the syntax: hspice -i Resistor_1.sp -o R_out
and look for the data in the R_out.lis file. For plotting the data refer step 2 in Section I above. Note that Resistor1.va and Resistor_1.sp should be \\
in the same folder.


Section III: Regarding MOSFETs

In NMOS the current flows from drain to source, so in the Verilog-A file use the format: I(d,s) <+ (expression for current) 
In PMOS the current flows from source to drain, so in the Verilog-A file use the format: I(s,d) <+ (expression for current) 
  





