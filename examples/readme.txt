
 RUN EXAMPLES (linux syntax):
   If MPC is installed in ~/Desktop/MPC/, 
   from the command line, with working directory of ~/Desktop/MPC/examples, 
   use the following command:

   1.)  'java -jar ../lib/MPC.jar example.mpc'

   This generates a JSim '.mod' ('example.mod') file that can added to a JSim project file:
   2.) 'jsim -f example.mod', then save new project file as needed.
   
   NOTE: currently, Working dir must be set to loction of .mpc file.

 Examples referenced in MPC15.10.05.Draft.pdf:
 - example.mpc
 - Ado2Ino3comp.mpc   <-- generate two substrate, mult-compartment model using ShortCodeLibrary.mod  


 HodgkinHuxley - Example showing the process of adding/removing channels and pumps as needed to accurately model sodium and potassium passage through cellular membrane.          
 
   If MPC is installed in ~/Desktop/MPC/, 
   from the command line, with working directory of ~/Desktop/MPC/examples/HodgkinHuxley, 
   use the following command:

   1.)  'java -jar ../../lib/MPC.jar hux_module.mpc'

   This generates a JSim '.mod' ('hux_module.mod') file that can added to a JSim project file:
   2.) 'jsim -f hux_module.mod', then save new project file as needed. 



