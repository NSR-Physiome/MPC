/* COPYRIGHT AND REQUEST FOR ACKNOWLEDGMENT OF USE:   
  Copyright (C) 2015 University of Washington. From the National Simulation Resource,  
  Director J. B. Bassingthwaighte, Department of Bioengineering, University of Washington, Seattle WA 98195-5061. 
  Academic use is unrestricted. Software may be copied so long as this copyright notice is included.

  This software was developed with support from NIH grants HL088516 and HL073598, NIBIB grant BE08417 
  and the Virtual Physiological Rat program GM094503 (PI: D.A.Beard). Please cite this grant in any 
  publication for which this software is used and send an email with the citation and, if possible, 
  a PDF file of the paper to: staff@physiome.org.  */

package ModConstruct;

import java.io.*;
import java.util.*;
class MPC {
    public static void main(String args[]) {
	try{
	    int num=args.length;
	    if( (num==1) ||(num==2) ) {
	    } else {
		System.out.println("Needs at least 1 argument: input fileName.");
		System.exit(0);
	    }
	    
	    if(!args[0].substring(args[0].length()-4,args[0].length()).equals(".mpc")) {
		System.out.println("MPC presource file must end in .mpc");
		System.exit(0);
	    }
	    String gather = new String();
	    // Open the file that is the first  command line parameter
	    FileInputStream istream = new FileInputStream(args[0]);
	    // Get the object of DataInputStream
	    DataInputStream in = new DataInputStream(istream);
	    BufferedReader br = new BufferedReader(new InputStreamReader(in));
	    String strLine;
	    // Prepare the FileWriter
	    String outfile = new String();
            if(num==2) { outfile=args[1];
	    } else {
	        outfile=args[0].substring(0,args[0].length()-4)+".mod";
	    }
	    FileWriter ostream = new FileWriter(outfile);
	    BufferedWriter out = new BufferedWriter(ostream);
	    // Create an ArrayList to hold the input file
	    ArrayList<String>inlines = new ArrayList<String>();
	    //Read File Line By Line
	    while ((strLine = br.readLine()) != null) {
		// add input line to ArrayList
		inlines.add(strLine);
	    }
	    //Close the input stream
	    in.close();
	    inlines.add("// This MML file generated from "+args[0]+" using MPC.");
	    
	    
	    ArrayList<String>gotCode = new ArrayList<String>();
	    boolean doneYet=false;
	    while(!doneYet) {  // Continue to process
		// Compress //%REPLICATE directives to single lines
		CompressDir bline = new CompressDir();
		ArrayList<String>replicate = new ArrayList<String>();
		replicate = bline.CompressDir(inlines);

                // Set mpc variables
                SetGlobalVal svline = new SetGlobalVal();
		ArrayList<String>replicat2 = new ArrayList<String>();
                replicat2 = svline.SetGlobalVal(replicate);

		// Replace code multiple times as required
		Replace cline = new Replace();
		ArrayList<String>compGet = new ArrayList<String>();
		compGet = cline.Replace(replicat2);
		
		// Get code from other models
		
		Get fline = new Get();
		
		//ArrayList<String>gotCode = new ArrayList<String>();
		gotCode = fline.Get(compGet);
		
		//May have gotten more REPLACE AND GET directives
		doneYet=true;
		for (int i=0; i<gotCode.size(); i++) {
		    String rline = gotCode.get(i);
		    if(rline.contains("//%REPLACE") ||
		    rline.contains("//%GET") )doneYet=false;
		}
		if(!doneYet) {
		    inlines.clear();
		    inlines.addAll(0,gotCode);
		}
	    } // end while
	    
	    // Remove exactly duplicated code lines
	    RemoveDup gline = new RemoveDup();
	    ArrayList<String> remdups = new ArrayList<String>();
	    remdups = gline.RemoveDup(gotCode);
	    
	    // Collect code with same names if specified, e.g.  "Ap:t"
	    Collect hline = new Collect();
	    ArrayList<String>collectedCode = new ArrayList<String>();
	    collectedCode = hline.Collect(remdups);
	    
	    // Insert new codeblock names
	    for (int i=0; i<collectedCode.size(); i++) {
		String s = collectedCode.get(i);
		if(s.contains("//%INSERTSTART")) collectedCode.set(i,s.replace("//%INSERTSTART","//%START"));
		if(s.contains("//%INSERTEND")) collectedCode.set(i,s.replace("//%INSERTEND","//%END"));
	    }
	    
	    // Write out assembled code
	    for (int i=0; i<collectedCode.size(); i++) {
		out.write( collectedCode.get(i) + "\n");
	    }
	    out.close();
	    
	} catch (Exception e){//Catch exception if any
	    System.err.println("Error: " + e.getMessage());
	}
    }
}
