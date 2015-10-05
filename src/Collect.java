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

// //%COLLECT("Name") Combine lines of code containing Name =
public class Collect  {
    public ArrayList<String>ilines;
    // constructor
    public ArrayList<String>Collect (ArrayList<String>ilines)  throws Exception {
	String lparen="(";
	String rparen=")";
	String quote="\"";
	String equals="=";
	String semicolon=";";
	String remove ="zZzzZZ";
	String oneComment="//";
	ArrayList<String>collect=new ArrayList<String>(ilines);
	
	// Find the Collect directive
	for (int i=0; i<collect.size(); i++) {
	    String s = collect.get(i);
	    // Locate //%COLLECT "item" directives
	    if (s.contains("//%COLLECT")) {
		String t = s;
		int ilptot = cntTimes(t,lparen);
		int irptot = cntTimes(t,rparen);
		int iqutot = cntTimes(t,quote);
		String msg = new String();
		if(ilptot!=irptot) {
		    msg = t + "  Unbalanced parentheses.  ";
		    throw new Exception(msg);
		}
		if (iqutot%2 !=0) {
		    msg = t + "  Odd number of quotes.  ";
		    throw new Exception(msg);
		}
		int nitems=iqutot/2;
		int ilp = s.indexOf(lparen);
		int istart = ilp;
		for (int items=0; items<nitems; items++) {
		    int is  = s.indexOf(quote,istart+1);
		    int ie  = s.indexOf(quote,is+1);
		    istart = ie;
		    // match is the item being collected
		    
		    String match = s.substring(is+1,ie).trim();
		    String offset = new String();
		    for (int ioff=0; ioff<match.length(); ioff++) {
			offset=offset+" ";
		    }
		    
		    // flag the collect directive with the remove string
		    if(items==0) {collect.set(i,remove);}
		    int jsave=-1;
		    int npieces=0;
		    ArrayList<String>code = new ArrayList<String>();
		    
		    // Scan all lines seeking equals sign
		    for (int j=0; j<collect.size(); j++) {
			String col = collect.get(j).trim();
			int ieq = col.indexOf(equals);
			
			if ( ieq >= 0 && (!col.contains(remove))) {
			    // Is the left hand side of the equation an exact match
			    // with Collect item?
			    if(match.equals(col.substring(0,ieq).trim())) {
				npieces++;
				if(npieces==1) {  // first collected piece
				    jsave=j;
				    String s1 =match+" = "+col.substring(ieq+1,col.length()).trim();
				    collect.set(j,remove);
				    code.add(s1);
				    while(!s1.contains(semicolon)) {
					npieces++;
					s1 = offset+" "+(String)collect.get(++j).trim();
					code.add(s1);
					collect.set(j,remove);
				    }
				} else { // additional collected pieces
				    
				    collect.set(j,remove);
				    col = col.substring(ieq+1,col.length()).trim();
				    if(!col.substring(0,1).equals("-") &&
				    !col.substring(0,1).equals("+") ) col = "+"+col;
				    col = offset+"   "+col;
				    code.add(col);
				    while(!col.contains(semicolon)) {
					npieces++;
					col =offset+(String)collect.get(++j).trim();
					collect.set(j,remove);
					
					code.add(col);
					
				    }
				}
			    }
			}
		    }
		    if(npieces>0 && jsave>=0) { // collect.set(jsave,composit);
			jsave--;
			for (int jcode=0; jcode<code.size(); jcode++) {
			    String f = code.get(jcode);
			    if (jcode<code.size()-1) {
				f =f.replaceAll(";","");
				collect.add(++jsave,f);
			    } else {
				if(!f.contains(";")) f=f+";";
				collect.add(++jsave,f);
			    }
			    
			}
			
		    }
		}
	    }
	}
	ArrayList<String>collected = new ArrayList<String>();
	for (int i=0; i<collect.size(); i++) {
	    if(!collect.get(i).contains(remove)) collected.add(collect.get(i));
	}
	return collected;
    }
    public static int cntTimes(String s1, String s2) {
	// Count number of times s2 in s1.
	int cnt = 0;
	int ndx = 0;
	while ((ndx = s1.indexOf(s2, ndx)) != -1) {
	    ++ndx;
	    ++cnt;
	}
	return cnt;
    }
    
    
}
