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

public class Get  {
	public ArrayList<String>ilines;
	// constructor
	public ArrayList<String>Get (ArrayList<String>ilines)  throws Exception {
	ArrayList<String>copy1 = new ArrayList<String>();
	int nlines = ilines.size();
	for (int i=0; i<nlines; i++) {
	    String s = (String)ilines.get(i);
	    if(s.contains("//%GET")) {
	    // Check for obvious errors
	    String lparen="(";
	    String rparen=")";
	    String quote="\"";
	    String equals="=";
	    int ilptot = cntTimes(s,lparen);
	    int irptot = cntTimes(s,rparen);
	    int iqutot = cntTimes(s,quote);
	    int ieqtot = cntTimes(s,equals);
	    String msg = new String();
	    if(ilptot!=irptot) {
	        msg = s + "  Unbalanced parentheses.  ";
	        throw new Exception(msg);
	    }
	    if (iqutot%2 !=0) {
	        msg = s + "  Odd number of quotes.  ";
	        throw new Exception(msg);
	    }
	    if(2*ieqtot-iqutot !=0) {
	        msg = s +"  Missing an equals sign or Extra quotes.";
	        throw new Exception(msg);
	    }
	    
	    ArrayList<String>extract = new ArrayList<String>();
	    extract = getCode(s);
	    for (int ic=0; ic<extract.size(); ic++) {
	        copy1.add((String)extract.get(ic));
	    }
	    } else {
	    copy1.add((String)s);
	    }
	}
	return copy1;
	}

	public static ArrayList<String> getCode (String g) {
	// PROCESS THE GET DIRECTIVE
	ArrayList<String>entity = new ArrayList<String>();
	ArrayList<String>code = new ArrayList<String>();
	int glength = g.length();
	// First substring is the fileName
	// Second substring is the codeBlockName
	// Then in pairs, "OldName=NewName" for
	// parameter/variable name replacements
	try {
	    int enddir = g.indexOf("//%GET")+6;
	    for ( int i=enddir; i<glength; i++) {
	    Character c = g.charAt(i);
	    int j = i;
	    if(entity.size()==0) {
	        if( isAlpha(c) || isLinux(c) ) {
	        while(!valSep(g.charAt(++i))) {}
	        entity.add(g.substring(j,i));
	        }
	    } else if(entity.size()==1) {
	        if(isAlpha(c)) {
	        while(!valSep(g.charAt(++i))) {}
	        entity.add(g.substring(j,i));
	        }
	    } else if( entity.size()%2==0) {
	        if(isQuote(c)) {
	        while(!isEquals(g.charAt(++i))) {}
	        entity.add(g.substring(j+1,i));
	        }
	    } else {
	        while(!isQuote(g.charAt(++i))) {}
	        entity.add(g.substring(j,i));
	    }
	    } // End separating strings
	    // Get file
	    if( entity.size()%2 !=0 ) {throw new Exception
	    ("Missing either fileName, codeBlockName, or part of replacement in //%GET."); }
	    String fileName = entity.get(0);
	    String codeBlockName = entity.get(1);
	    boolean atStart = false;
	    boolean atEnd=false;
	    // Open the file that is the first  command line parameter
	    FileInputStream istream = new FileInputStream(fileName);	    
	    // Get the object of DataInputStream
	    DataInputStream in = new DataInputStream(istream);
	    BufferedReader br = new BufferedReader(new InputStreamReader(in));
	    String s;
	    //Read File Line By Line
	    while ((s = br.readLine()) != null) {
	    // add input line to ArrayList

            // Need to differentiate between //%START SodiumChannel and //%START ExternalSodiumChannel
		String[] words = s.split("\\s+"); // split s into separate identifiers based on whitespace
		if ( s.contains("//%START") && (words.length >1) && words[1].equals(codeBlockName) )  {
			atStart=true;
			continue;
	        }
		if ( s.contains("//%END") && (words.length >1) && words[1].equals(codeBlockName) ) {
		    atEnd = true;
		}
	    
	    if(atStart && !atEnd) {
	        // PROCESS CODE; ADD TO OUTPUT
	        int nentities=(entity.size()-2)/2;
	        for (int i=0; i<nentities; i++) {
	        int ibegin=0;
		//Get old name and new name
	        String oldNam= entity.get(2*i+2);
	        String newNam= entity.get(2*i+3);
		//remove any leading or trailing whitespace
		oldNam = oldNam.trim();
		newNam = newNam.trim();
	        while(  s.indexOf(oldNam, ibegin) > -1) {
	            int j = s.indexOf(oldNam, ibegin);
	            int iend = s.length();
	            int io = oldNam.length();
	            int inm = newNam.length();
	            String r = new String();
	            // Is the oldNam bounded by either valid
	            // separators or starts or ends the line of code
	            if( (j==0) && (io==iend) ) {//s is oldnam
	            r=newNam;
	            ibegin = inm-1;
	            s = r;
	            if(inm==1) {ibegin+=1;}
	            }
	            else if( (j==0) && (io<iend) && validSeparator(s.charAt(j+io))) 
			{// s is oldname +more
			    r=newNam+s.substring(j+io,iend);
			    ibegin = j+io;
			    s=r;
			}
	            else if(  (j>0) && (j+io<iend) && validSeparator(s.charAt(j-1)) && validSeparator(s.charAt(j+io)))
			{ // s is more+oldNam+more
			    r = s.substring(0,j)+newNam+s.substring(j+io,iend);
			    ibegin = j+inm-1;
			    if(inm==1) {ibegin+=1;}
			    s = r;
			}
	            else if ( (j>0) && (j+io==iend) && validSeparator(s.charAt(j-1)) ) 
			{ //s is more+oldNam
			    r=s.substring(0,j)+newNam;
			    ibegin = j+inm;
			    s = r;
			}
	            else { ibegin++;}
	        }
	        }
	        if(atStart && !s.contains("//%START") || !s.contains(codeBlockName)) {
		    code.add(s);
	        }
	    }
	    }
	    String err = "In file "+fileName+" missing either //%START or //%END with codeBlockName="+codeBlockName+" ";
	    if( !atStart || !atEnd ) throw new Exception (err);
	    
	    //Close the input stream
	    in.close();
	} catch (Exception e){//Catch exception if any
	    System.err.println("Error: " + e.getMessage());
	}
	return code;
	}
	
	public static boolean valSep( Character c) {
	String separator=" ()=,";
	String tab = "\t";
	Character TAB = tab.charAt(0);
	if(TAB.equals(c) ) return true;
	for (int i=0; i<separator.length(); i++) {
	    Character sep = separator.charAt(i);
	    if (sep.equals(c) ) return true;
	}
	return false;
	}
	
	public static boolean isLinux( Character c ) {
	String separator=".~*?[]-!{}";
	for (int i=0; i<separator.length(); i++) {
	    Character sep = separator.charAt(i);
	    if (sep.equals(c) ) return true;
	}
	return false;
	}
	
	
	public static boolean validSeparator( Character c) {
	String tab = "\t";
	Character TAB = tab.charAt(0);
	if(TAB.equals(c) ) return true;
	String separator=" {}[]()^*/+-=><@,;:.\"";
	for (int i=0; i<separator.length(); i++) {
	    Character sep = separator.charAt(i);
	    if (sep.equals(c) ) return true;
	}
	return false;
	}
	
	public static boolean isAlpha (Character c) {
	Character a = new Character('a');
	Character z = new Character('z');
	Character A = new Character('A');
	Character Z = new Character('Z');
	if ((c>=a && c<=z) || (c>=A && c<=Z)) {
	    return true;
	} else {
	    return false;
	}
	
	}
	
	public static boolean isAlpha_ ( Character c ) {
	if(isAlpha(c) || c.equals('_') ) {
	    return true;
	} else {
	    return false;
	}
	}
	
	public static boolean isNum (Character c ) {
	Character zero = new Character('0');
	Character nine = new Character('9');
	if (c>=zero && c<=nine) {
	    return true;
	} else {
	    return false;
	}
	}
	
	public static boolean isAlpha_Num ( Character c ) {
	if (isAlpha(c) || isNum(c) ) {
	    return true;
	} else {
	    return false;
	}
	}
	
	public static boolean isQuote( Character c) {
	if(c.equals('"'))  {
	    return true;
	} else {
	    return false;
	}
	}
	public static boolean isEquals( Character c) {
	if(c.equals('='))  {
	    return true;
	} else {
	    return false;
	}
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



