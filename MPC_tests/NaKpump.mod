// MODEL NAME: Winslow1999_INaK_module
// SHORT DESCRIPTION: This Sodium-Potassium pump comes from Winslow etal. 1999 model that analyzes 
// the influence of voltage-dependent Ca2+ independent transient current (Ito1) on the action potential 
// duration in normal vs failing canine and human cardiac myocytes. 


import nsrunit;
unit conversion on;

math   INaK {


// NaK pump from Winslow Rice Jafri 1999, flip currents for hodgkin-Hux consistency
//%START timeDomain
  realDomain t msec; t.min=0; t.max=100; t.delta=0.10; // time
//%END timeDomain

//%START universalConsts
    real   F    = 96.5 coulomb/mmol,           
           Temp = 279.5 K,       // 6.3 C in HH squid experiments
	   R    = 8.314 J/(mol*K),
           Cm  = 1 uF/cm^2;
    real   zK =1 dimensionless, zNa =1 dimensionless;   // Valency of K, Na ion, unitless   
//%END universalConsts

//%START NaKpump
//%START externalInputsSodiumPotPump
// Comment out/delete as needed:
// For HodgkinHuxley model the inward current is positive.
    extern real minusI(t) uA/cm^2;
    extern real Vclamp(t) mV;            
    real clamp_0no_1yes = 0; 
    real Vdepolar  = -95.0 mV;
    real VV(t) mV, 
          V(t) mV; // displacement of the membrane potential from rest (depolarization negative);
    when (t=t.min) { VV  = if (clamp_0no_1yes=0) Vdepolar  else Vclamp; }
    V   = if (clamp_0no_1yes = 0) VV else Vclamp;
    real Kinit_o  = 10.0 mM,                     
	 Nainit_o = 600 mM,
	 Kinit_i  = 15.0 mM,                    
	 Nainit_i = 10.0 mM; 
    real K_o(t) mM, K_i(t) mM;
    real Na_o(t) mM, Na_i(t) mM;

    real Vcell = 3.8E-5 uL;            // Cell volume
    real Vfcyto = 0.68 dimensionless;  // cyto vol fraction: cyto vol per cell vol (ml/ml)
    real Vcyto uL;                     // Cytoplasm volume
         Vcyto = Vcell * Vfcyto; 
//%END externalInputsSodiumPotPump

  real     Acap = 1.534*10^(-4) cm^2,       // membrane capacitance area
	   Vmyo = 25.84*10^(-6) uL,
           NaKmax = 0.693 uA/uF,            // Max pump velocity
           INaKmax uA/cm^2;
           INaKmax = Cm*NaKmax;
  real     KmNa_i = 10.0 mM,
           KmK_o  = 1.5 mM;     
 
  real     INaK(t) uA/cm^2, fNaK(t), sigma(t);
  real     INapump(t) uA/cm^2; 
  real     IKpump(t)  uA/cm^2; 

  INaK   = INaKmax*fNaK/(1+(KmNa_i/Na_i)^1.5)     // The output is INaK.
               * K_o/(K_o+KmK_o);                            //(A.51)
  fNaK   = 1/(1+0.1245*exp(-0.1*(V*F/(R*Temp)))
               +0.0365*sigma*exp(-(V*F/(R*Temp))));          //(A.52)
  sigma  = (exp(Na_o/(67.3 mM))-1)/7;                  

  INapump   = -3*INaK; // flip sign from orig for HodgkinH modeling
  IKpump    = 2*INaK;  // flip sign from orig for HodgkinH modeling

//%END NaKpump

//  NEEDED to make module a stand alone model.
 real Iion(t) uA/cm^2;  // Total ionic current
      Iion    = INaK;
      VV:t   = if (clamp_0no_1yes = 0) (minusI - Iion)/Cm   else  0;  

 when (t=t.min) { K_i = Kinit_i; Na_i = Nainit_i; }
     K_o = Kinit_o;    // Assume extracellular conc is const 
     Na_o = Nainit_o;  // Assume extracellular conc is const 

 Na_i:t = -(Acap/(zNa*Vcyto*F))*( INapump);
 K_i:t = -(Acap/(zK*Vcyto*F))*(IKpump);
//  END of stand alone section

} //end

/*

 DETAILED DESCRIPTION:
 Sodium-Potassium pump module from Winslow et al. 1999 paper. Assumes extracellular 
 Na+ and K+ concentrations are constant.
 
 SHORTCOMINGS/GENERAL COMMENTS:
	- This is a stand-alone NaK pump taken from Winslow et al. 1999 paper.
        - Current sign flipped to make compatible with Hodgkin-Huxley squid model.
 
 KEY WORDS: cardiac action potential, NaK pump, Sodium-Potassium pump, MPC Module, 
 Cell Physiology, PMID10082479   

 REFERENCES:
 Winslow R.L., Rice J., Jafri S., Marban E., and O'Rourke B., Mechanisms of altered 
 excitation-contraction coupling in canine tachycardia-induced heart failure, II: model 
 studies. Circ Res. 1999 Mar 19;84(5):571-86.

 Greenstein J.L., Wu R., Po S., Tomaselli G.F., and Winslow, R.L., Role of the Calcium-independent 
 transient outward current Ito1 in shaping action potential morphology and duration. 
 Circ Res, 1026-1033, 2000.


 REVISION HISTORY:
	Original Author : JBB  Date: 10/sep/15
        Revised by: BEJ Date:30sep15 : Make into NaKpump module

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



