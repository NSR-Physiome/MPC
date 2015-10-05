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
// //%SETGLOBALVAL %name%= ("newname")
//   Relace each occurence of %name% with newname
//   replicating the line in which %name% occurs.
//   This is usually used to set Numbers in Replace directives,
//   e.g.
//   //%SETGLOBALVAL %N%=("5")
//   //%REPLACE (%im1%=("#0#%N%-2"), %i%=("#1#%N%-1"), %ip1%=("#2#%N%") )
//   c0:t = 0;
//   c%i%:t = D*(c%im1%-2*c%i% +c%ip1%)/dx^2;
//   c%N% =0;
//   //%ENDREPLACE
//   produces
//   c0:t = 0;
//   c1:t = D*(c0-2*c1 +c2)/dx^2;
//   c2:t = D*(c1-2*c2 +c3)/dx^2;
//   c3:t = D*(c2-2*c3 +c4)/dx^2;
//   c4:t = D*(c3-2*c4 +c5)/dx^2;
//   c5:t = 0;
//   The scope of of SetGlobalVar is the entire .mpc file.
//   SetGlobalVal is executed after the compress and before the replace directives
public class SetGlobalVal  {
    public ArrayList<String>ilines;
    // constructor
    public ArrayList<String>SetGlobalVal (ArrayList<String>ilines)  throws Exception {
	ArrayList<String> copy1 = new ArrayList <String>(ilines);
	boolean done=false;
	while (!done) {
	    int nlines = copy1.size();
	    for (int i=0; i<nlines; i++) {
		String t = (String)copy1.get(i);
		if( t.contains("//%SETGLOBALVAL")) {
// Check for obvious errors
		    String lparen="(";
		    String rparen=")";
		    String quote="\"";
		    String equals="=";
		    String percent="%";
		    int ilptot = cntTimes(t,lparen);
		    int irptot = cntTimes(t,rparen);
		    int iqutot = cntTimes(t,quote);
		    int ieqtot = cntTimes(t,equals);
		    int iprtot = cntTimes(t,percent);
		    String msg = new String();
		    if(ilptot!=irptot) {
			msg = t + "  Unbalanced parentheses.  ";
			throw new Exception(msg);
		    }
		    if (iqutot%2 !=0) {
			msg = t + "  Odd number of quotes.  ";
			throw new Exception(msg);
		    }
		    if(iprtot%2 !=1) {
			// includes the % in //%REPLACE
			msg = t +"  Missing a % sign.";
			throw new Exception(msg);
		    }
		    ArrayList<String> copy2 = new ArrayList<String>(nlines);
		    int n0pct = t.indexOf(percent);
		    int n1pct = t.indexOf(percent,n0pct+1);
		    int n2pct = t.indexOf(percent,n1pct+1);
		    String sname = t.substring(n1pct,++n2pct).trim();
		    
		    int n1quot = t.indexOf(quote);
		    int n2quot = t.indexOf(quote,n1quot+1);
		    String sname1 = t.substring(++n1quot,n2quot).trim();
		    
		    for (int j=0; j<i; j++) {
			String s = (String)copy1.get(j);
			if (s.contains(sname)) { s=s.replaceAll(sname,sname1);
			}
			copy2.add(j,s);
		    }
  		    String s = "  ";
		    copy2.add(i,s);
		    for (int j=i+1; j<nlines; j++) {
			s = (String)copy1.get(j);
			if (s.contains(sname)) { s=s.replaceAll(sname,sname1);
			}
                        int k=j - 1;
			copy2.add(k,s);
		    }
		    copy1.clear();
		    copy1.addAll(copy2);
		}
	    }
	    done=true;
	}
	return copy1;
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

