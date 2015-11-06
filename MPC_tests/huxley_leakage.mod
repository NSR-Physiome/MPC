// MODEL NUMBER: 
// MODEL NAME: HodgkinHuxley_LeakageCurrent
// SHORT DESCRIPTION: Hodgkin and Huxley (HH 1952d): Nerve action potential for squid giant 
// axon. Quantitative module of time and voltage-dependent transmembrane leak current, Ileak. 



// Hodgkin-Huxley Model

import  nsrunit;
unit conversion on;

math HodgkinHuxleyleakage {

//%START timeDomain
  realDomain t msec; t.min=0; t.max=100; t.delta=0.10; // time
//%END timeDomain

//%START universalConsts
 real Cm       = 1 uF/cm^2;
 real Temp     = 279.5 K;      // 6.3 Celsius in HH squid experiments
 real Vnorm    = 1 mV; // for non-dimensionalizing V in function expressions, 
//%END universalConsts

//%START LeakageCurrent
//%START externalInputsLeakCurrent
// For HH model the inward current is positive.
    extern real minusI(t) uA/cm^2;
    extern real Vclamp(t) mV;            
    real clamp_0no_1yes = 0; 
    real Vdepolar = -90 mV;   
    real VV(t) mV, 
          V(t) mV;   // displacement of the membrane potential from rest (depolarization negative);
    when (t=t.min) { VV = if (clamp_0no_1yes=0) Vdepolar  else Vclamp; }
    V   = if (clamp_0no_1yes = 0) VV else Vclamp;
//%END externalInputsLeakCurrent

// *** Leakage current ***
 real gbar0    = 0.3 mmho/cm^2,  // max leak conductance
      Vl       = -10.613 mV;
 real Il(t) uA/cm^2;  // Ionic currents
 Il      = gbar0 * (V-Vl);
//%END LeakageCurrent

// Needed for stand alone model:
  real Iion(t) uA/cm^2;
       Iion  = Il;  
  VV:t = if (clamp_0no_1yes = 0) (-minusI-Iion)/Cm else 0; // universal var
}

/*

 DETAILED DESCRIPTION:
 The Hodgkin-Huxley model of membrane current is concerned with the flow of electric currents through 
 the surface membrane of a nerve fiber. The total current is given by the first equation (see Equations) 
 and consists of the capacity current,Cm, dV/dt, and an ionic current, split into three component 
 currents carried by sodium ions (INa), potassium ions (IK) and other ions (Il). The ionic 
 currents are written in terms of conductances, gNa, gK, and gbarl.	

 SHORTCOMINGS/GENERAL COMMENTS:
	- Specific inadequacies or next level steps
 
 KEY WORDS: squid nerve action potential, ionic currents, sodium, potassium, voltage clamp, 
 neuron, Publication, Cell physiology, PMID2185861, MPC module 

 REFERENCES:
  Hodgkin AL and Huxley AF. A quantitative description of membrane current and its application 
  to conduction and excitation in nerve. J Physiol. 500-544, 1952d.

  Cole KS and Moore JW. Potassium ion current in the squid giant axon: dynamic characteristic. 
  Biophys J 1: 1-14, 1960.

  Hille B. Ionic Channels of Excitable Membranes. Sunderland, MA: Sinauer, 2001.

  Jack JJB, Noble D, and Tsien RW. Electric Current Flow in Excitable Cells. Oxford, England: 
  Clarendon Press, 1975.

 REVISION HISTORY:
	Original Author : JBB  Date: 06/12/01
	Revised by: BEJ Date:14nov11 : Update comment format
	Revised by: BEJ Date:28jun12 : Add keyword reproducible
        Revised by: BEJ Date:31jul13 : Update comments, formatting

 COPYRIGHT AND REQUEST FOR ACKNOWLEDGMENT OF USE:   
  Copyright (C) 1999-2015 University of Washington. From the National Simulation Resource,  
  Director J. B. Bassingthwaighte, Department of Bioengineering, University of Washington, Seattle WA 98195-5061. 
  Academic use is unrestricted. Software may be copied so long as this copyright notice is included.
  
  This software was developed with support from NIH grant HL073598. 
  Please cite this grant in any publication for which this software is used and send an email 
  with the citation and, if possible, a PDF file of the paper to: staff@physiome.org. 

*/

