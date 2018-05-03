/* SHORT DESCRIPTION: A modified version of A2B with
   relabeling the two species for adenosine (Ado)
   and inosine (Ino) in three regions (plasma (p), 
   interstitial fluid region (i) and parenchymal cell 
   (c) ) with addition of competitive transporter for 
   Ado and Ino on cell membrane and enzyme conversion 
   of Ado to Ino.
*/
import nsrunit; unit conversion on;
math Ado2Ino { 
//%START ado2inoModel
// INDEPENDENT VARIABLES
realDomain  t s; t.min=0; t.max=30; t.delta = 0.1;
real L = 0.1 cm;
real Ndivx = 31;
realDomain x cm ; x.min=0; x.max=L; x.ct = Ndivx;

// PARAMETERS
real Flow = 1 ml/(g*min); // Flow rate 
real PSgAdo = 3 ml/(g*min); // Exchg rate 
real PSpcAdo = 5 ml/(g*min); // Exchg rate 
real PSgIno = 3 ml/(g*min); // Exchg rate 
real PSpcIno = 5 ml/(g*min); // Exchg rate 
real Gado2ino = 10 ml/(g*min); // Conversion rate 
real Vp = 0.05 ml/g; // Volume of Vp
real Vi = 0.15 ml/g; // Volume of Vi
real Vc = 0.60 ml/g; // Volume of Vc
real DpAdo = 1e-6 cm^2/sec; // Diffusion coeff 
real DiAdo = 1e-6 cm^2/sec; // Diffusion coeff 
real DcAdo = 1e-6 cm^2/sec; // Diffusion coeff 
real DpIno = 1e-6 cm^2/sec; // Diffusion coeff 
real DiIno = 1e-6 cm^2/sec; // Diffusion coeff 
real DcIno = 1e-6 cm^2/sec; // Diffusion coeff 
extern real Adoin(t) mM; // Inflowing concentration
extern real Inoin(t) mM; // Inflowing concentration
// DEPENDENT VARIABLES
real Adoout(t) mM; // Outflowing concentration 
real Inoout(t) mM; // Outflowing concentration 
real Adop(t,x) mM; // Concentration 
real Adoi(t,x) mM; // Concentration 
real Adoc(t,x) mM; // Concentration 
real Inop(t,x) mM; // Concentration 
real Inoi(t,x) mM; // Concentration 
real Inoc(t,x) mM; // Concentration 
// INITIAL CONDITIONS
when(t=t.min) Adop=0;
when(t=t.min) Adoi=0;
when(t=t.min) Adoc=0;
when(t=t.min) Inop=0;
when(t=t.min) Inoi=0;
when(t=t.min) Inoc=0;
// BOUNDARY CONDITIONS
when  (x=x.min) (-Flow*L/Vp)*(Adop-Adoin)+DpAdo*Adop:x = 0; 
when  (x=x.max) { Adop:x = 0; Adoout = Adop;}
when  (x=x.min) (-Flow*L/Vp)*(Inop-Inoin)+DpIno*Inop:x = 0; 
when  (x=x.max) { Inop:x = 0; Inoout = Inop;}
when(x=x.min)  Adoi:x=0;
when(x=x.max)  Adoi:x=0; 
when(x=x.min)  Adoc:x=0;
when(x=x.max)  Adoc:x=0; 
when(x=x.min)  Inoi:x=0;
when(x=x.max)  Inoi:x=0; 
when(x=x.min)  Inoc:x=0;
when(x=x.max)  Inoc:x=0; 
// ADDITIONAL PARAMETERS
private real ZEROT = 0 mmol/cm^2;// Removes Initial Conditions from Parameter List
private real ZEROM = 0 mM;	// Removes Initial Conditions from Parameter List
//     ENZYME CONVERSION PARAMETERS
real Etot = 0.1 mM;
real kf1 = 2 mM^(-1)*s^(-1);
real kb1 = 1 s^(-1);
real kf2 = 1 s^(-1);
real kb2 = 1 mM^(-1)*s^(-1);
//     COMPETITIVE TRANSPORTER PARAMETERS
real Ttot	= 7e-6 mmol/cm^2;	// Transporter density on membrane
real KdAdoi	= 1 mM;		// Equilib Dissoc const 
real KdAdoc	= 1 mM;		// Equilib Dissoc const 
real KdInoi	= 1 mM;		// Equilib Dissoc const 
real KdInoc	= 1 mM;		// Equilib Dissoc const 
real konAdoi	= 10 mM^(-1)*s^(-1);	// Bind rate 
real konAdoc	= 10 mM^(-1)*s^(-1);	// Bind rate 
real konInoi	= 10 mM^(-1)*s^(-1);	// Bind rate 
real konInoc	= 10 mM^(-1)*s^(-1);	// Bind rate 
real kofAdoi	= KdAdoi*konAdoi;		// Dissoc rate 
real kofAdoc	= KdAdoc*konAdoc;		// Dissoc rate 
real kofInoi	= KdInoi*konInoi;		// Dissoc rate 
real kofInoc	= KdInoc*konInoc;		// Dissoc rate 
real S 		= 1 cm^2/g;		// Surface area
real SoVi 	= S/Vi;		// Surface to volume ratio
real SoVc 	= S/Vc;		// Surface to volume ratio
real kTi2c 	= 100 sec^(-1);	// Flip rate 
real kTc2i 	= 100 sec^(-1);	// Flip rate 
real kTAdoi2c	 = 100 sec^(-1);// Flip rate ;
real kTAdoc2i	 = 100 sec^(-1);// Flip rate ;
real kTInoi2c	 = 100 sec^(-1);// Flip rate ;
real kTInoc2i	 = 100 sec^(-1);// Flip rate ;
real Enz(t,x) mM;		// Concentration of Enz in Vc
real ECmplx(t,x) mM;		// Concentration of ECmplx in Vc
real TAdoi(t,x) mmol/cm^2;		// Transporter complex
real TAdoc(t,x) mmol/cm^2;		// Transporter complex
real TInoi(t,x) mmol/cm^2;		// Transporter complex
real TInoc(t,x) mmol/cm^2;		// Transporter complex
real Ti(t,x) mmol/cm^2;		// Free transporter
real Tc(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz	= Etot;
when(t=t.min) ECmplx	= ZEROM;
when(t=t.min) TAdoi 	= ZEROT;
when(t=t.min) TAdoc 	= ZEROT;
when(t=t.min) TInoi 	= ZEROT;
when(t=t.min) TInoc 	= ZEROT;
when(t=t.min) Ti	= Ttot/2;
when(t=t.min) Tc	= Ttot/2;
when(x=x.min)  Enz:x=0;
when(x=x.max)  Enz:x=0; 
when(x=x.min)  ECmplx:x=0;
when(x=x.max)  ECmplx:x=0; 
when(x=x.min)  TAdoi:x=0;
when(x=x.max)  TAdoi:x=0; 
when(x=x.min)  TAdoc:x=0;
when(x=x.max)  TAdoc:x=0; 
when(x=x.min)  TInoi:x=0;
when(x=x.max)  TInoi:x=0; 
when(x=x.min)  TInoc:x=0;
when(x=x.max)  TInoc:x=0; 
when(x=x.min)  Ti:x=0;
when(x=x.max)  Ti:x=0; 
when(x=x.min)  Tc:x=0;
when(x=x.max)  Tc:x=0; 
// PDE CALCULATIONS
Adop:t = -(Flow*L/Vp)*Adop:x
     + DpAdo*Adop:x:x 
       +PSgAdo/Vp*(Adoi-Adop);
Inop:t = -(Flow*L/Vp)*Inop:x
     + DpIno*Inop:x:x 
       +PSgIno/Vp*(Inoi-Inop);
Adoi:t = DiAdo*Adoi:x:x
       +PSgAdo/Vi*(Adop-Adoi)
       +PSpcAdo/Vc*(Adoc-Adoi)
       +(-konAdoi*Adoi*Ti + kofAdoi*TAdoi)*SoVi;
Inoi:t = DiIno*Inoi:x:x
       +PSgIno/Vi*(Inop-Inoi)
       +PSpcIno/Vc*(Inoc-Inoi)
       +(-konInoi*Inoi*Ti + kofInoi*TInoi)*SoVi;
TAdoi:t = konAdoi*Adoi*Ti - kofAdoi*TAdoi
        - kTAdoi2c*TAdoi + kTAdoc2i*TAdoc;
TAdoc:t = konAdoc*Adoc*Tc - kofAdoc*TAdoc
        +kTAdoi2c*TAdoi - kTAdoc2i*TAdoc;
TInoi:t = konInoi*Inoi*Ti - kofInoi*TInoi
        - kTInoi2c*TInoi + kTInoc2i*TInoc;
TInoc:t = konInoc*Inoc*Tc - kofInoc*TInoc
        +kTInoi2c*TInoi - kTInoc2i*TInoc;
Ti:t = -konAdoi*Adoi*Ti + kofAdoi*TAdoi
     -konInoi*Inoi*Ti + kofInoi*TInoi
     - kTi2c*Ti + kTc2i*Tc;
Tc:t = -konAdoc*Adoc*Tc + kofAdoc*TAdoc
     -konInoc*Inoc*Tc + kofInoc*TInoc
     +kTi2c*Ti - kTc2i*Tc;
Adoc:t = DcAdo*Adoc:x:x
       +PSpcAdo/Vc*(Adoi-Adoc)
       -Gado2ino/Vc*Adoc
       +(-konAdoc*Adoc*Tc + kofAdoc*TAdoc)*SoVc
       -kf1*Adoc*Enz + kb1*ECmplx;
Inoc:t = DcIno*Inoc:x:x
       +PSpcIno/Vc*(Inoi-Inoc)
       +Gado2ino/Vc*Adoc
       +(-konInoc*Inoc*Tc + kofInoc*TInoc)*SoVc
       +kf2*ECmplx - kb2*Inoc*Enz;
Enz:t = -(kf1*Adoc + kb2*Inoc)*Enz
      + (kb1 + kf2)*ECmplx;
ECmplx:t = (kf1*Adoc + kb2*Inoc)*Enz
         - (kb1 + kf2)*ECmplx;
//%END ado2inoModel
}
// This MML file generated from Ado2Ino.mpc using MPC v1.01.
