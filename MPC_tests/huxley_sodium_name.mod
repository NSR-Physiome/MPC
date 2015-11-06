// MODEL NUMBER: 
// MODEL NAME: HodgkinHuxley_SodiumChannel
// SHORT DESCRIPTION: Hodgkin and Huxley (HH 1952d): Nerve action potential for squid giant 
// axon. Sodium channel part of whole model is here. Time and voltage-dependent transmembrane 
//  currents for Na+. 



// Hodgkin-Huxley Na+ channel 

import  nsrunit;
unit conversion on;


math HodgkinHuxleySodium {

//%START timeDomain
  realDomain t msec; t.min=0; t.max=100; t.delta=0.10; // time
//%END timeDomain

//%START SodiumChanneluniversalConsts
 real Cm       = 1 uF/cm^2;
 real Temp     = 279.5 K,      // 6.3 Celsius in HH squid experiments
      F        = 96.5 coulomb/mmol,  // Faraday const         
      R        = 8.314 J/(mol*K),
      phi      = 3^((Temp-279.5)/(283- 273 K)); // Orig, in Celsius: ((Temp-6.3)/10)
 real Vnorm    = 1 mV; // for non-dimensionalizing V in function expressions, 
//%END SodiumChanneluniversalConsts

//%START externalInputsSodiumChannel
// For HH model the inward current is positive.
    extern real minusI(t) uA/cm^2;
    extern real Vclamp(t) mV;            
    real clamp_0no_1yes = 0;
    real Vdepolar = -90 mV; 
    real VV(t) mV, 
          V(t) mV;   // displacement of the membrane potential from rest (depolarization negative);
    when (t=t.min) { VV = if (clamp_0no_1yes=0) Vdepolar  else Vclamp;  }
    V   = if (clamp_0no_1yes = 0) VV else Vclamp;   

    real Nainit_o = 600 mM;   // Init extramembrane conc, seawater
    real Nainit_i = 7 mM;     // Init innermembrane conc: No calc off of init nernst
    real Na_o(t) mM, Na_i(t) mM;
    real VNa(t) mV;           // Nernst potential
    VNa = -(R*Temp/F)*ln(Na_o/Na_i);   // Nernst eq
//%END externalInputsSodiumChannel

//%START SodiumChannel
// *** Sodium channel *** 
  real gbarNa   = 120 mmho/cm^2;  // max sodium conductance 
  
  real alpham0 = 2.5/(exp(2.5)-1)*(1  msec^(-1)),
       betam0  = 4 msec^(-1),
       alphah0 = 0.07 msec^(-1),
       betah0  = 1/(exp(3)+1) *(1 msec^(-1));

     //Variables for all the algebraic equations
  real INa(t) uA/cm^2,  // Ionic currents
       alpham(t) msec^(-1),   // rate constant of activating molecules from out to in
       betam(t) msec^(-1),    // rate constant of activating molecules from in to out
       alphah(t) msec^(-1),   // rate constant of inactivating molecules from out to in
       betah(t) msec^(-1),    // rate constant of inactivating molecules from in to out
       gNa(t) mmho/cm^2;  // Sodium conductance

  real minf(t), hinf(t);
  private real taum(t) msec, tauh(t) msec;
  // State variables for all the ODEs
  real  m(t),   // proportion of activating molecules on the inside
        h(t);   // proportion of inactivating molecules on the outside

   when (t=t.min) {
       m = alpham0/(alpham0+betam0);
       h = alphah0/(alphah0+betah0);
     }

   alpham  = (1 msec^-1)*(if (V/Vnorm = -25) 1
              else 0.1*(V/Vnorm+25)/(exp((V/Vnorm+25)/10)-1));
   betam   = (1 msec^-1)*(4*exp((V/Vnorm)/18));
   alphah  = (1 msec^-1)*(0.07*exp((V/Vnorm)/20));
   betah   = (1 msec^-1)*(1/(exp(( V/Vnorm+30)/10)+1));
   minf    = alpham/(alpham+betam);
   hinf    = alphah/(alphah+betah);
   tauh    =1/(alphah+betah);
   taum    =1/(alpham+betam);
   gNa     = gbarNa * m^3 * h;
   INa     = gNa * (V-VNa);
   m:t = phi*(alpham*(1-m)-betam*m);
   h:t = phi*(alphah*(1-h)-betah*h);

//%END SodiumChannel

// Needed for stand alone model:
  real Iion(t) uA/cm^2;
       Iion = INa;
  VV:t = if (clamp_0no_1yes = 0) (-minusI-Iion)/Cm else 0; // universal var
  Na_i = Nainit_i;  // assume const
  Na_o = Nainit_o;

}

/*

 DETAILED DESCRIPTION:
 The Hodgkin-Huxley model of membrane current is concerned with the flow of electric currents through 
 the surface membrane of a nerve fiber. The total current is given by the first equation (see Equations) 
 and consists of the capacity current,Cm, dV/dt, and an ionic current, split into three component 
 currents carried by sodium ions (INa), potassium ions (IK) and other ions (Il). The ionic 
 currents are written in terms of conductances, gNa, gK, and gbarl.	

 SHORTCOMINGS/GENERAL COMMENTS:
   - This is a module to be added to larger models as needed. It can be run as a stand-alone
     model as well.
 
 KEY WORDS: squid nerve action potential, ionic currents, sodium, potassium, voltage clamp, 
 Nernst, neuron, Cell physiology, PMID2185861, MPC module 

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
        Revised by: BEJ Date:30sep15 : Make into Na+ channel module

 COPYRIGHT AND REQUEST FOR ACKNOWLEDGMENT OF USE:   
  Copyright (C) 1999-2015 University of Washington. From the National Simulation Resource,  
  Director J. B. Bassingthwaighte, Department of Bioengineering, University of Washington, Seattle WA 98195-5061. 
  Academic use is unrestricted. Software may be copied so long as this copyright notice is included.
  
  When citing JSim please use this reference: Butterworth E, Jardine BE, Raymond GM, Neal ML, Bassingthwaighte JB. 
  JSim, an open-source modeling system for data analysis [v3; ref status: indexed, http://f1000r.es/3n0] 
  F1000Research 2014, 2:288 (doi: 10.12688/f1000research.2-288.v3)  

  This software was developed with support from NIH grants HL088516 and HL073598, NIBIB grant BE08417 
  and the Virtual Physiological Rat program GM094503 (PI: D.A.Beard). Please cite this grant in any 
  publication for which this software is used and send an email with the citation and, if possible, 
  a PDF file of the paper to: staff@physiome.org. 

*/

