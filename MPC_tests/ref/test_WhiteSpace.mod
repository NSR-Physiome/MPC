// Build up Hinch 2004 model with the following modules:
// Tab in mod file  block name line causes error



import nsrunit;
unit conversion on;
unit mJ_per_mole_K=.001 kilogram^1*meter^2*second^(-2)*kelvin^(-1)*mole^(-1);
unit ms=.001 second^1;
unit um3_per_ms=1E-15 meter^3*second^(-1);
unit per_ms=1E3 second^(-1);
unit per_ms3=1E9 second^(-3);
unit mM=1 meter^(-3)*mole^1;
unit um3=1E-18 meter^3;
unit C_per_mole=1 second^1*ampere^1*mole^(-1);
unit um3_mM_per_ms=1E-15 second^(-1)*mole^1;
unit mM_per_ms=1E3 meter^(-3)*second^(-1)*mole^1;

math hinch2004modules {
	realDomain t ms; t.min=0; extern t.max; extern t.delta;
    real Cm = 1 uF/cm^2;  // membrane cap per area
    real F C_per_mole;    // Faraday constant
    F=96487;
    real F_Ca = 2*F;      // Faraday const, Ca2+ 
    real R mJ_per_mole_K; // Universal gas constant
    R=8314.5; 
    real T kelvin;
    T=295;
    extern real Ca_o mM;   // Extracellular calcium (2+) concentration
    extern real Ca_i mM;   // Intracellular calcium (2+) concentration
    extern real Ca_SR mM;  // Calcium (2+) concentration in sarcoplasmic reticulum
    extern real minusI(t) uA/cm^2;
    extern real Vclamp(t) mV;            
    real clamp_0no_1yes = 0; 
    real Vdepolar = -90 mV;   
    real VV(t) mV, 
          V(t) mV;   // displacement of the cell membrane potential from rest (depolarization negative);
    when (t=t.min) { VV = if (clamp_0no_1yes=0) Vdepolar  else Vclamp; }
    V   = if (clamp_0no_1yes = 0) VV else Vclamp;
    real FVRT(t) dimensionless;    // Convenience variable that is equal to F*V/(R*T)
    real FVRT_Ca(t) dimensionless; // Convenience variable that is equal to (F*V/(R*T))*2
 

// Put the currents together and external conc dependencies, etc:
// Ca2+ current going into the cytoplasm is positive.
  real Iion(t) uA/cm^2;
       Iion  = ((I_RyR+I_LCC+ I_TRPN+ I_SR - I_pCa+ I_CaB -I_SERCA)*F_Ca* - (I_NaCa*F))*V_myo/Area_mem; // I_NaCa: 2*F?
  VV:t = if (clamp_0no_1yes = 0) (-minusI-Iion)/Cm else 0; // universal var
Ca_i:t = betaCDMA*(I_RyR+I_LCC+ I_TRPN+ I_SR -I_NaCa - I_pCa+ I_CaB -I_SERCA );
Ca_SR:t = (V_myo/VSR)* (-RyR +I_SERCA - I_SR); // signs??
}
// This MML file generated from test_WhiteSpace.mpc using MPC v1.0.
