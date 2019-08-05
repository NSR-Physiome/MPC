// Test directive names can be found if one contains the other
import  nsrunit;
unit conversion on;

math huxley_modules {
  realDomain t msec; t.min=0; t.max=100; t.delta=0.10; // time
 real Cm       = 1 uF/cm^2;
 real Temp     = 279.5 K,      // 6.3 Celsius in HH squid experiments
      F        = 96.5 coulomb/mmol,  // Faraday const         
      R        = 8.314 J/(mol*K),
      phi      = 3^((Temp-279.5)/(283- 273 K)); // Orig, in Celsius: ((Temp-6.3)/10)
 real Vnorm    = 1 mV; // for non-dimensionalizing V in function expressions, 
//%START externalInputsKChannel
// For HH model the inward current is positive.
    extern real minusI(t) uA/cm^2;
    extern real Vclamp(t) mV;            
    real clamp_0no_1yes = 0; 
    real Vdepolar = -90 mV; 
    real VV(t) mV, 
          V(t) mV;   // displacement of the membrane potential from rest (depolarization negative);
    when (t=t.min) { VV = if (clamp_0no_1yes=0) Vdepolar  else Vclamp; }
    V   = if (clamp_0no_1yes = 0) VV else Vclamp;  

    real Kinit_o = 10 mM;       // Init extramembrane conc, seawater
    real Kinit_i = 150 mM;      // Init innermembrane conc: No calc off of init nernst
    real K_o(t) mM, K_i(t) mM;  // K values outer and inner
    real VK(t) mV;                  // Nernst potential
    VK = -(R*Temp/F)*ln(K_o/K_i);   // Nernst eq
//%END externalInputsKChannel

//
// *** Potassium channel *** 
     real   gbarK    = 36 mmho/cm^2;   // max potassium conductance      
     real Kalphan0 = 0.1/(exp(1)-1)*(1 msec^(-1)),  // always use V=0
          Kbetan0  = 0.125 msec^(-1);               // to calculate i.c.

     // Variables for all the algebraic equations
     real IK(t) uA/cm^2,  // Ionic currents
          Kalphan(t) msec^(-1),   // rate constant of particles from out to in
          Kbetan(t) msec^(-1),    // rate constant from in to out 
          gK(t) mmho/cm^2;       // potassium conductance
     real Kninf(t);
     real Ktaun(t) msec;
     
     real Kn(t);   // proportion of the particles in a certain position
      // IC:
      when (t=t.min) { Kn = Kalphan0/(Kalphan0+Kbetan0); }

      Kalphan  = (1 msec^-1)*(if (V/Vnorm = -10) 0.1
                  else   0.01*(V/Vnorm+10)/(exp((V/Vnorm+10)/10)-1));
      Kbetan   = (1 msec^-1)*(0.125*exp((V/Vnorm)/80));
      Kninf    = Kalphan/(Kalphan+Kbetan);
      Ktaun    =1/(Kalphan+Kbetan);
      gK      = gbarK  * Kn^4;
      IK      = gK  * (V-VK);
      Kn:t = phi*(Kalphan*(1-Kn)-Kbetan*Kn);

    real clamp_0no_1yes = 0;
    when (t=t.min) { VV = if (clamp_0no_1yes=0) Vdepolar  else Vclamp;  }
    V   = if (clamp_0no_1yes = 0) VV else Vclamp;   

    real Nainit_o = 600 mM;   // Init extramembrane conc, seawater
    real Nainit_i = 7 mM;     // Init innermembrane conc: No calc off of init nernst
    real Na_o(t) mM, Na_i(t) mM;
    real VNa(t) mV;           // Nernst potential
    VNa = -(R*Temp/F)*ln(Na_o/Na_i);   // Nernst eq
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



}




// This MML file generated from huxley_name_contains.mpc using MPC v1.02.
