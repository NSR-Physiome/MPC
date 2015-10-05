
package ModConstruct;

import java.io.*;
import java.util.*;
import java.lang.*;
import javax.script.ScriptEngineManager;
import javax.script.ScriptEngine;
// //%REPLACE %name%= ("name1","name2" ...)
//   relace each occurence of %name% with name1, name 2 ...
//   replicating the line in which %name% occurs.
//   RESTRICTIONS: name1, name2, ... cannot contain "="
public class Replace  {
    public ArrayList<String>lines;
    // constructor
    public ArrayList<String>Replace (ArrayList<String>ilines)  throws Exception {
	ArrayList <String> copy1 = new ArrayList <String>();
	ArrayList <String> copy3 = new ArrayList <String>();
// Supply missing ENDREPLACE directives at end of file
	int nlines = ilines.size();
        int cntRandER=0;
        for (int i=0; i<nlines;i++) {
            copy3.add((String)ilines.get(i));
            if( ((String)copy3.get(i)).contains("//%REPLACE")) {
               cntRandER+=1;}
            if( ((String)copy3.get(i)).contains("//%ENDREPLACE")) {
               cntRandER+=-1; }
        }
        if(cntRandER>0) {
            String ER="//%ENDREPLACE";
            for (int i=0; i<cntRandER; i++) {
                copy3.add(ER);
            }
        }
        nlines=copy3.size();
	int nDoCnt=0;
	int nEndDoCnt=0;
	int maxEmbedCnt=0;
	int embedCnt=0;
// Determine the levels of embeddedness
	for (int i=0; i<nlines; i++) {
	    copy1.add((String)copy3.get(i));
	    if( ((String)copy1.get(i)).contains("//%REPLACE")) {
		String t = (String)copy1.get(i);
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
		nDoCnt++;
		embedCnt++;
		maxEmbedCnt=Math.max(embedCnt,maxEmbedCnt);
	    }
	    if( ((String)copy1.get(i)).contains("//%ENDREPLACE")) {
		nEndDoCnt++;
		embedCnt--;
		maxEmbedCnt=Math.max(embedCnt,maxEmbedCnt);
	    }
	}
	// Count number of embedded constructs at each level
	int recount = -1;
	int numEmbed[] = new int[maxEmbedCnt];
	for (int i=0; i<nlines; i++) {
	    if( ((String)copy1.get(i)).contains("//%REPLACE")) {
		recount+=1;
		numEmbed[recount]++;
	    }
	    if( ((String)copy3.get(i)).contains("//%ENDREPLACE")) {
		recount-=1;
	    }
	}
	
	if(nDoCnt!=nEndDoCnt) {
	    String msg="# of REPLACEs, "+nDoCnt+", not equal # of ENDREPLACEs, "+nEndDoCnt+".";
	    throw new Exception(msg);
	}
	ArrayList <String> copy2  = new ArrayList <String>();
	ArrayList <String> hold   = new ArrayList <String>();
	ArrayList <String> expand = new ArrayList <String>();
// Process the most embedded REPLACE directives first	
// Note the structure used, maxEmbedCnt is decremented inside the loop
	for (int jEmbedCnt=0; jEmbedCnt<maxEmbedCnt; ) {
	    int ncopy1=copy1.size();
	    int traceReplaceEndR =0;
	    for (int i=0; i<ncopy1; i++) {
		String s = (String)copy1.get(i);
		if (s.contains("//%REPLACE")) {traceReplaceEndR +=1;  }
		if (s.contains("//%ENDREPLACE")) {traceReplaceEndR += -1;  }
		if (s.contains("//%REPLACE")) {
		    if (traceReplaceEndR==maxEmbedCnt)   {
			hold.clear();
			hold.add(s);
			while(!((String)copy1.get(++i)).contains("//%ENDREPLACE")) {
			    hold.add((String)copy1.get(i));
			}
			
			expand = replicate(hold);
			
			int nexpand = expand.size();
			for (int iexpand=0; iexpand<nexpand; iexpand++) {
			    copy2.add((String)expand.get(iexpand));
			}
			expand.clear();
			traceReplaceEndR += -1;
		    } else {
			copy2.add(s);
		    }
		} else {
		    copy2.add(s);
		}
	    }
	    copy1.clear();
	    int ncopy2=copy2.size();
	    for (int i=0; i<ncopy2; i++) {
		copy1.add((String)copy2.get(i));
	    }
	    copy2.clear();
	    numEmbed[maxEmbedCnt-1]+=(-1);
	    if(numEmbed[maxEmbedCnt-1]<=0) --maxEmbedCnt;
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
    
    public static ArrayList <String> replicate (ArrayList <String> hold) throws Exception {
	int nhold = hold.size();
	String s = ExpandN((String)hold.get(0));
	// Decode REPLACE directive
	int nvar = cntTimes(s,"=");
	int nrepl = cntTimes(s,"\"")/(2*nvar);
	s=s.replaceAll("//%REPLACE","     ");
	String[] var = new String[nvar];
	String[][] repl = new String[nvar][nrepl];
	for (int ivar=0; ivar<nvar; ivar++) {
	    // Get replacement names
	    var[ivar]=getvar(s);
	    s=s.replaceFirst(var[ivar]," ");
	    s=s.replaceFirst("="," ");
	    
	    String[] grepl = new String[nrepl];
	    // Get replacement values
	    grepl = getrepl(s, nrepl);
	    for (int irepl=0; irepl<nrepl; irepl++) {
		repl[ivar][irepl] = grepl[irepl];
		s=s.replaceFirst(repl[ivar][irepl]," ");
		s=s.replaceFirst("\""," ");
		s=s.replaceFirst("\""," ");
	    }
	}
	// Allow for maximum string replacements
	String[] newhold = new String[nrepl*(nhold-1)];
	int scnt=0;
	for (int i=1; i<nhold; i++) { //1st string was REPLACE directive after compression
	    String hs = (String)hold.get(i);
	    boolean replacement=false;
	    for (int ivar=0; ivar<nvar; ivar++) {
		if (hs.contains(var[ivar])) replacement=true;
	    }
	    if(replacement) {
		for (int irepl=0; irepl<nrepl; irepl++) {
		    newhold[scnt] = hs;
		    for (int ivar=0; ivar<nvar; ivar++) {
			newhold[scnt]=newhold[scnt].replaceAll(var[ivar],repl[ivar][irepl]);
		    }
		    scnt+=1;
		}
	    } else {
		newhold[scnt] = hs;
		scnt+=1;
	    }
	}
	
	ArrayList<String>expanded = new ArrayList<String>();
	for (int i=0; i<scnt; i++) {
	    expanded.add(newhold[i]);
	}
	return expanded;
    }
    
    public static String ExpandN (String s) throws Exception {
	// //%REPLACE name1%=("#2#12")  Replace name1 with 2,3,4,...11, and 12
	String quote="\"";
	String Nsym="#";
	int expansion = cntTimes(s,Nsym);
	if(expansion>0) {
	    int two = 2;
	    int times = expansion/two;
	    for (int i=0; i<times; i++) {
		int firstN=s.indexOf(Nsym);
		int secondN=s.indexOf(Nsym,firstN+1);
		int rightquote=s.indexOf(quote,secondN+1)+1;
		int leftquote= firstN;
		while(!s.substring(leftquote,leftquote+1).equals(quote)) {leftquote--;}
		ScriptEngineManager mgr = new ScriptEngineManager();
		ScriptEngine engine = mgr.getEngineByName("JavaScript");
		
		String ms1 = s.substring(firstN+1,secondN);
		String s1 = engine.eval(ms1).toString();
		double ds1 = Double.valueOf(s1);
		int i1 = (int)ds1;
		
		String ms2 = s.substring(secondN+1,rightquote-1);
		String s2 = engine.eval(ms2).toString();
		double ds2 = Double.valueOf(s2);
		int i2 = (int)ds2;
		
		String name=s.substring(leftquote+1,firstN).trim();
		String augment = new String();
		if(i1<=i2) {
		    for ( int it=i1; it<i2+1; it++) {
			augment =augment+" \""+name+it+"\", ";
		    }
		} else {
		    for ( int it=i1; it>i2-1; it--) {
			augment =augment+" \""+name+it+"\", ";
			
		    }
		}
		int iend = s.length();
		String sbegin = s.substring(1,leftquote);
		String send   = s.substring(rightquote,iend);
		s=sbegin+augment+send;
		
	    }
	}
	return s;
    }
    
    
    public static int back2pct(String s, int eqls) {
	Character pc = new Character('%');
	int pccnt = 0;
	for (int i=eqls; i>0; i--) {
	    Character is = s.charAt(i);
	    if ( (is.equals(pc)) ){
		pccnt++;
		if (pccnt==2) {
		    return i;
		}
	    }
	}
	return -1;
    }
    
    
    public static String getvar (String s) {
	String var = new String();
	int eqls  = s.indexOf("=");
	int start = back2pct(s,eqls);
	if(start>=0 && start<=eqls) {
	    var =  s.substring(start, eqls).trim() ;
	}
	return var;
    }
    
    public static String[] getrepl (String s, int nrepl ) {
	// locate quoted strings
	String s1 = new String(s);
	String[] repl = new String[nrepl];
	for (int irepl=0; irepl<nrepl; irepl++) {
	    int start = s1.indexOf("\"");
	    s1=s1.replaceFirst("\""," ");
	    int end   = s1.indexOf("\"");
	    repl[irepl]=s1.substring(start,end).trim();
	    s1=s1.replaceFirst("\""," ");
	}
	return repl;
    }
}

