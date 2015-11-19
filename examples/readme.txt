
* TO RUN MPC EXAMPLES (linux syntax):
   If MPC is installed in ~/Desktop/MPC/, 
   from the command line, with working directory of ~/Desktop/MPC/examples, 
   use the following command:

   1.)  'java -jar ../lib/MPC.jar example_intro.mpc'

   This generates a JSim '.mod' ('example_intro.mod') file that can added to a JSim project file:
   2.) 'jsim -f example.mod', then save new project file as needed.
   
   NOTE: currently, Working dir must be set to location of .mpc file.


 *******************************************
 Example referenced in F1000 software paper (Submitted Nov 2015):
 - F1000_example.mpc <-- Example.mpc from paper.
   - F1000_CodeLibrary.mod <--- Library used by F1000_example.mpc to generate JSim model.

 Generate JSim model:
   1. 'java -jar ../lib/MPC.jar F1000_example.mpc' --> generates 'F1000_example.mod' file
   2. Then open up with JSim (http://www.physiome.org/jsim/): 'jsim -f F1000_example.mod'

 - F1000_example.proj  <--- Two compartment MPC example imported and made to run in JSim with plot (fig 2 of F1000 paper).
   To run, install JSim (http://www.physiome.org/jsim/), then open F1000_example.proj with JSim. 
  
 *******************************************
 Examples referenced in MPC_details_date.pdf:
 - example_intro.mpc  <- Intro MPC code to generate simple JSim model.
 - Ado2Ino3comp.mpc   <-- generate two substrate, mult-compartment JSim model using ShortCodeLibrary.mod  
   - ShortCodeLibrary.mod <-- Code library used by Ado2Ino3comp.mpc
   - A3comp.mod  <-- Model used by Ado2Ino3comp.mpc

 *******************************************
 HodgkinHuxley - Example showing the process of adding/removing channels and pumps as needed to accurately model sodium and potassium passage through cellular membrane.          
 
   If MPC is installed in ~/Desktop/MPC/, 
   from the command line, with working directory of ~/Desktop/MPC/examples/HodgkinHuxley, 
   use the following command:

   1.)  'java -jar ../../lib/MPC.jar hux_module.mpc'

   This generates a JSim '.mod' ('hux_module.mod') file that can added to a JSim project file:
   2.) 'jsim -f hux_module.mod', then save new project file as needed. 



