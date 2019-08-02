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
import java.lang.*;
// from the beginning of a 
// math ModelName {
// section to the end,
// locate exactly duplicated statements and delete them
// retaining statements that contain comment delimiters
//   
public class RemoveDup  {
    public ArrayList<String>ilines;
    
    // constructor
    public ArrayList<String>RemoveDup (ArrayList<String>ilines) throws Exception {
	ArrayList<String>copy1 = new ArrayList<String>(ilines);
	ArrayList<String>copy2 = new ArrayList<String>();
	String remove="zZzzZZZZ";
	boolean mathflag=false;
        int deleteBlockCnt = 0;
        for (int i=0; i<ilines.size()-1; i++) {
            if (deleteBlockCnt<0) {
                String msg = new String();
                msg = "  Misplaced //%ENDDELETE directive.  ";
                    throw new Exception(msg);
                }
            String d1 = ilines.get(i);
            if(d1.contains("//%STARTDELETE")) {
                deleteBlockCnt +=1;
            }
            if(deleteBlockCnt>0) {copy1.set(i,remove); }
            else {copy1.set(i,ilines.get(i));}
            if(d1.contains("//%ENDDELETE")) {
                deleteBlockCnt -= 1; 
                copy1.set(i,remove); 
            }
        }
        if (deleteBlockCnt!=0) {
                String msg = new String();
                msg = "  # of //%ENDDELETE directives not equal to # of //%STARTDELETE directives  ";
                    throw new Exception(msg);
                }

	for (int i=0; i<ilines.size()-1; i++) {
	    String s1=ilines.get(i);
	    String s1trim = s1.trim();
            if( s1trim.startsWith("DIAGRAM") ) {mathflag=false;}
            if( s1trim.startsWith("DETAILED")) {mathflag=false;}

	    if( (!mathflag) && (s1trim.startsWith("math ") ) ) {
		mathflag=true;
	    } else {
		if (mathflag) {
		    if(s1trim.contains("//%COM") )
		    {copy1.set(i,remove); }
		    else if(!s1trim.equals("/*") &&
		    !s1trim.equals("*/") &&
		    !s1trim.equals("//") )  {
			boolean notdone = true;
			for (int j=i+1; j<ilines.size(); j++) {
			    String s2= ilines.get(j);
			    String s2trim = s2.trim();
                            if( (mathflag) && (notdone) && (s2trim.startsWith("DIAGRAM") ) ) {notdone=false;}
                            if( (mathflag) && (notdone) && (s2trim.startsWith("DETAILED") ) ) {notdone=false;}
			    if(s1.equals(s2) && mathflag && notdone && !s1trim.startsWith("+") && !s1trim.startsWith("-" )) { 
 
                            if (s2.trim().length()>0) { // Do not delete blank lines
								copy1.set(j,remove);  
								System.out.println("deleting "+s2+" ");}
                            }
			}
		    }
		}
	    }
	}
	for (int i=0; i<ilines.size(); i++) {
            if( copy1.get(i).contains("//%COM")) {copy1.set(i,remove);}
	    if(!copy1.get(i).contains(remove) ) copy2.add(copy1.get(i)) ;
	}
	return copy2;
    }
}
