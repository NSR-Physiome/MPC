/* COPYRIGHT AND REQUEST FOR ACKNOWLEDGMENT OF USE:   
  Copyright (C) 2015 University of Washington. From the National Simulation Resource,  
  Director J. B. Bassingthwaighte, Department of Bioengineering, 
  University of Washington, Seattle WA 98195-5061. 
  Academic use is unrestricted. Software may be copied so long as this copyright notice is included.

  This software was developed with support from NIH grants HL088516 and HL073598, NIBIB grant BE08417 
  and the Virtual Physiological Rat program GM094503 (PI: D.A.Beard). Please cite this grant in any 
  publication for which this software is used and send an email with the citation and, if possible, 
  a PDF file of the paper to: staff@physiome.org.  */

package ModConstruct;

import java.io.*;
import java.util.*;
import java.lang.*;
// Assemble continued directive on a single line
public class CompressDir  {
    public ArrayList<String>ilines;
    // constructor
    public ArrayList<String> CompressDir (ArrayList<String> ilines)  throws Exception {
	ArrayList <String> copy1 = new ArrayList <String>();
	int nlines = ilines.size();
	String dir1 = "//%REPLACE";
	String dir2 = "//%GET";
	String dir3 = "//%COLLECT";
	boolean anotherDir=false;
	// Compress continued directives into single line
	for (int i=0; i<nlines; i++) {
	    String s = (String)ilines.get(i);
	    if (s.contains(dir1) || s.contains(dir2) || s.contains(dir3))  { // This is a directive
               anotherDir=false;

		while (parensMatchCnt(s)  ) {
		    String sCont =(String)ilines.get(++i);
		    int iCont=sCont.indexOf("//%");

                    if( sCont.indexOf("//%START")          >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%END")            >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%INSERTSTART")    >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%INSERTEND")      >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%SETGLOBALVAL")   >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%RELACE")         >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%ENDREPLACE")     >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%GET")            >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%COLLECT")        >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%COM")            >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%STARTDELETE")    >=0) {anotherDir=true;}
                    if( sCont.indexOf("//%ENDDELETE")      >=0) {anotherDir=true;}
                    if(anotherDir) {
			throw new Exception
			("Missing continuation directive \n" +s);
		    } else {
			s=s.trim()+" "+sCont.replaceAll("//%","   ").trim();
		    }
		}
		copy1.add(s);
	    } else {
		copy1.add(s);
	    }
	}
	return copy1;
    }
    
    public static boolean parensMatchCnt (String s) {
	// Left and Right Parentheses must be equal and there must
	// be at least one left parenthesis
	int leftp =cntTimes(s,"(");
	int rightp =cntTimes(s,")");
	int currentCnt = +leftp-rightp;
        boolean done=false;
        boolean doContinue=true;
        if( (currentCnt==0) && (leftp==0)) {return doContinue;}
	if( (currentCnt==0) && (leftp>0) ) {return done;}
        if(  currentCnt>0                ) {return doContinue;}
        if(  currentCnt<0                ) {return done;}
	return false;
    }
    
    public static int cntTimes(String s1, String s2) {
	// Count number of times s2 in s1.
	// But do not count \s2 (s2 preceded by back slash
	 /*   String backslash="\\";
	    Character BACKSLASH=backslash.CharAt(0); */
	int cnt = 0;
	int ndx = 0;
	while ((ndx = s1.indexOf(s2, ndx)) != -1) {
	    ++ndx;
	    ++cnt;
	}
	return cnt;
    }

  }


