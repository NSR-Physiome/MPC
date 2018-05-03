/* MODEL NAME MultiflowAdo2Ino
   SHORT DESCRIPTION: 3 region model for Adenosine and
   inosine with competitive transporter on cell membrane
   and enzyme conversion of adenosine to inosine in the
   cell in a 10 path heterogeneous flow model.
*/
//*****************************************************************************
/* Statistics on a curve */

source procedure curveStatistics(input@t,output@t; ai,ti,RDi,ao,tao,RDo,tsys,RDsys) {

/* ai  (ao)       area of input (output) curve
   ti  (to)       transit time of input (output) curve
   RDi (RDo)      relative dispersion of input (output) curve
   tsys           transit time of the system
   RDsys          relative dispersion of the system

   NOTA BENE, the results are not dimensionally scaled or checked.
   A basic assumption is that both the input and output have the
   same units.  

       language="java";
       maincode={{
           RegularGridData t = (RegularGridData) input.grid(0);
           double S1=0; double S2=0; double S3=0;
           double T1=0; double T2=0; double T3=0;
           for (int i=0; i<t.ct(); i++) {
               double time = t.min()+((double)i)*t.delta();
               double time2 = time*time;
               double timeO=time;
               double timeO2=timeO*timeO;
               double in=input.realVal(i);
               double ot=output.realVal(i);
               S1+=in;
               S2+=in*time;
               S3+=in*time2;
               T1+=ot;
               T2+=ot*timeO;
               T3+=ot*timeO2;
           }
           ai.set(t.delta()*S1);
           double tI=S2/S1;
           ti.set(tI);
           double rdI=( Math.sqrt((Math.abs(S3/S2)/(S2/S1))-1));
           RDi.set(rdI);

           ao.set(t.delta()*T1);
           double tO=T2/T1;
           tao.set(tO);
           double rdO=( Math.sqrt((Math.abs(T3/T2)/(T2/T1))-1));
           RDo.set(rdO);

           tsys.set(tO-tI);
           RDsys.set( Math.sqrt(rdO*tO*rdO*tO-rdI*tI*rdI*tI)/(tO-tI) );
         }};
} // END
//
//
/*
   DETAILED DESCRIPTION
// 
// ALGORITHM FOR HETEROGENEOUS FLOWS
   1.) An algorithm for heterogeneous flow:
   A probability density function (PDF) is generated. Suggested
   choices are LagNormal or Gaussian, with area = 1, tMean =1,
   and RD = .4 to .55. The PDF represents the distribution of
   the input concentration among the various paths.
//
   2.) A range of relative flows is selected, relFmin, and relFmax.
   Reasonable values are relFmin=0.2; relFmax=2.0. The user
   chooses the number of flow paths, PATHS. Reasonable values
   are 5 to 20. The reader chooses how the relative flows are
   generated: (RECOMMENDED) equally spaced in the logarithm of the 
   relative flows; equally spaced in the transit time of the relative flows
   (emphasizes the slowest speeds);  or equally spaced in the relative flows
   (preserves the shape of the PDF but ignores the slowest flows). 

   ADVICE ON CHOOSING THE RELATIVE FLOW MINIMUM
   fmin should be less than the value in this chart when using the
   Lagged Normal function with tMean=1
                      skewness    
      rd              0             .5          1        1.5         2       
      .3             .5          .4683      .4324      .3912     .3437   
      .35            .5          .4626      .4193      .3684     .3076   
      .4             .5          .4568      .4056      .3437     .2674   
      .45            .5          .4509      .3912      .3170     .2221   
      .5             .5          .4449      .3762      .2881     .1710   
      .55            .5          .4387      .3603      .2566     .1126   
      .6             .5          .4324      .3437      .2221     .0454 
//
   3.) The OUTPUT for the user is wts(NPATH), and normf(NPATH).
//
   4.) This program is designed to function with the MODULAR PROGRAM
   CONSTRUCTOR. The .mpc file using this code must be modified and 
   run everytime the number of relative flow paths is changed. 
//
   5.) Individual path flows should be summed as
       Cout(i)*wts(i)*normf(i)
   6.) To use this routine change the number of flows, i.e.
       "10=7" for 7 flows
   7.) To get the weighted flow in each path
ChemOutWeightedPathNumber = ChemOutPathNumber*wts(PathNumber)*normf(PathNumber)
ChemOutTot = ChemOutWeightedPathnumber

MPC's COLLECT(ChemoutTot)  will collect the total outflow

*/
import nsrunit; unit conversion on;
math MultiflowAdo2Ino{

real PATHS = 10; // SET PATHS in presource file
//
   private real NPATH   = PATHS*1;
   public real NPATHS   = NPATH;
   real relFmin = 0.2 ;
   real relFmax = 2.0 ;
//
   realDomain NP ; NP.min=1; NP.max=NPATH; NP.delta=1;
   private NP.min, NP.max, NP.delta, NP.ct;
   private real userF(NP), userWt(NP);
// THESE ARE THE USER DEFINED FLOWS AND WEIGHTS
// SEE CHOICES UNDER SPACING
real UserF1   = 0;
real UserF2   = 0;
real UserF3   = 0;
real UserF4   = 0;
real UserF5   = 0;
real UserF6   = 0;
real UserF7   = 0;
real UserF8   = 0;
real UserF9   = 0;
real UserF10   = 0;
real UserWt1  = 0;
real UserWt2  = 0;
real UserWt3  = 0;
real UserWt4  = 0;
real UserWt5  = 0;
real UserWt6  = 0;
real UserWt7  = 0;
real UserWt8  = 0;
real UserWt9  = 0;
real UserWt10  = 0;
userF =
        if( abs(NP-1)<0.1) UserF1 else
        if( abs(NP-2)<0.1) UserF2 else
        if( abs(NP-3)<0.1) UserF3 else
        if( abs(NP-4)<0.1) UserF4 else
        if( abs(NP-5)<0.1) UserF5 else
        if( abs(NP-6)<0.1) UserF6 else
        if( abs(NP-7)<0.1) UserF7 else
        if( abs(NP-8)<0.1) UserF8 else
        if( abs(NP-9)<0.1) UserF9 else
        if( abs(NP-10)<0.1) UserF10 else
        0;
userWt =
         if( abs(NP-1)<0.1) UserWt1 else
         if( abs(NP-2)<0.1) UserWt2 else
         if( abs(NP-3)<0.1) UserWt3 else
         if( abs(NP-4)<0.1) UserWt4 else
         if( abs(NP-5)<0.1) UserWt5 else
         if( abs(NP-6)<0.1) UserWt6 else
         if( abs(NP-7)<0.1) UserWt7 else
         if( abs(NP-8)<0.1) UserWt8 else
         if( abs(NP-9)<0.1) UserWt9 else
         if( abs(NP-10)<0.1) UserWt10 else
         0;
//
   realDomain NPp1 ; NPp1.min=1; NPp1.max=NPATH+1; NPp1.delta=1;
   private NPp1.min, NPp1.max, NPp1.delta, NPp1.ct;
//
   realDomain rflo ; rflo.min=relFmin; rflo.max=relFmax; rflo.delta=0.01;
   private rflo.min, rflo.max,  rflo.delta, rflo.ct;
// SET PDF FOR RELATIVE FLOWS
   extern real relFpdf(rflo) ;	// pdf for relative flows
// CALCULATE CUMULATIVE PDF
   real cumpdf(rflo) ;
   when(rflo=rflo.min) cumpdf=0;
   cumpdf:rflo=relFpdf;
// CHOOSE SPACING FOR RELATIVE FLOWS
   choice SPACING ("= in log(relative flow)","= in transit time","= in relative flow","User specifies") = 1;
   private real spacing = SPACING;
//
// SPACE FLOW BINS EVENLY IN SPACING DOMAIN
   private real srelFmin , srelFmax , swidth ;
   srelFmin = if(spacing=1) ln(relFmin) else
                if(spacing=2) 1/relFmin   else
                relFmin;
   srelFmax = if(spacing=1) ln(relFmax) else
                if(spacing=2) 1/relFmax   else
                relFmax;
   swidth = (srelFmax-srelFmin)/NPATH;
//
// COMPUTE ENDS OF EACH RELATIVE FLOW BIN
   real ends(NPp1) ;         // ends of each flow bin
   ends(NPp1) = if(spacing=1) exp(srelFmin+(NPp1-1)*swidth) else
             if(spacing=2) 1/ (srelFmin+(NPp1-1)*swidth) else
                           (srelFmin+(NPp1-1)*swidth);
//
// MAKE EACH RELATIVE FLOW THE CENTER OF ITS FLOW BIN
   real f(NP) ;       // relative flows
   f(NP) = if(spacing<>4) (ends(NP+1)+ends(NP))/2 else userF(NP);
//
// INTERPOLATE cumpdf(rflo) TO GET cumValues(NPp1);
   real cumValues(NPp1) ;    // cumpdf(rflo)->cumValues(PATHS+1)
   cumValues(NPp1) = cumpdf(ends(NPp1));
//
// GENERATE WT FOR EACH BIN so that integral wt*df=1
   real wt(NP) ;             // d(cumValues)/d(ends)
   wt(NP)=if(spacing<>4) (cumValues(NP+1)-cumValues(NP))/(ends(NP+1)-ends(NP)) else userWt(NP);
   real wdsum ;              // for normalization
   wdsum = sum( NP=NP.min to NP.max,wt);
   real wxfsum;
   wxfsum = sum( NP=NP.min to NP.max, wt*f);
   real flownormal = wdsum/wxfsum;
//
// NORMALIZE WEIGHTS
   real wts(NP) ;            // normalized weights
   wts(NP) = wt(NP)/wdsum;
   real normf(NP);
   normf = f*flownormal;
// SET MEAN FLOW
   real Fmean = 1 ml/(g*min);
   real Fi(NP) ml/(g*min);
   Fi=Fmean*normf;
   real Fi1 ml/(g*min);
   real Fi2 ml/(g*min);
   real Fi3 ml/(g*min);
   real Fi4 ml/(g*min);
   real Fi5 ml/(g*min);
   real Fi6 ml/(g*min);
   real Fi7 ml/(g*min);
   real Fi8 ml/(g*min);
   real Fi9 ml/(g*min);
   real Fi10 ml/(g*min);
// CALCULATE FLOW FOR EACH PATH
   Fi1=Fi(1);
   Fi2=Fi(2);
   Fi3=Fi(3);
   Fi4=Fi(4);
   Fi5=Fi(5);
   Fi6=Fi(6);
   Fi7=Fi(7);
   Fi8=Fi(8);
   Fi9=Fi(9);
   Fi10=Fi(10);
//% returns Fi1, Fi2 ... FiN
//
// INDEPENDENT VARIABLES
realDomain  t s; t.min=0; t.max=30; t.delta = 0.1;
real L = 0.1 cm;
real Ndivx = 31;
realDomain x cm ; x.min=0; x.max=L; x.ct = Ndivx;
//
// PARAMETERS
real PSgAdo = 3 ml/(g*min);	// Exchg rate 
real PSpcAdo = 5 ml/(g*min);	// Exchg rate 
real PSgIno = 3 ml/(g*min);	// Exchg rate 
real PSpcIno = 5 ml/(g*min);	// Exchg rate 
real Gado2ino = 10 ml/(g*min);	// Conversion rate 
real Vp = 0.05 ml/g;		// Volume of Vp
real Vi = 0.15 ml/g;		// Volume of Vi
real Vc = 0.60 ml/g;		// Volume of Vc
real DpAdo = 1e-6 cm^2/sec;	// Diffusion coeff 
real DiAdo = 1e-6 cm^2/sec;	// Diffusion coeff 
real DcAdo = 1e-6 cm^2/sec;	// Diffusion coeff 
real DpIno = 1e-6 cm^2/sec;	// Diffusion coeff 
real DiIno = 1e-6 cm^2/sec;	// Diffusion coeff 
real DcIno = 1e-6 cm^2/sec;	// Diffusion coeff 
extern real Adoin(t) mM;	// Inflowing concentration
extern real Inoin(t) mM;	// Inflowing concentration
// DEPENDENT VARIABLES
real Adoout1(t) mM;		// Outflowing concentration 
real Inoout1(t) mM;		// Outflowing concentration 
real Adop1(t,x) mM;		// Concentration 
real Adoi1(t,x) mM;		// Concentration 
real Adoc1(t,x) mM;		// Concentration 
real Inop1(t,x) mM;		// Concentration 
real Inoi1(t,x) mM;		// Concentration 
real Inoc1(t,x) mM;		// Concentration 
// INITIAL CONDITIONS
when(t=t.min)  Adop1=0;
when(t=t.min)  Adoi1=0;
when(t=t.min)  Adoc1=0;
when(t=t.min)  Inop1=0;
when(t=t.min)  Inoi1=0;
when(t=t.min)  Inoc1=0;
// BOUNDARY CONDITIONS
when  (x=x.min) (-Fi1*L/Vp)*(Adop1-Adoin)+DpAdo*Adop1:x = 0; 
when  (x=x.max) { Adop1:x = 0; Adoout1 = Adop1;}
when  (x=x.min) (-Fi1*L/Vp)*(Inop1-Inoin)+DpIno*Inop1:x = 0; 
when  (x=x.max) { Inop1:x = 0; Inoout1 = Inop1;}
when(x=x.min)  Adoi1:x=0;
when(x=x.max)  Adoi1:x=0; 
when(x=x.min)  Adoc1:x=0;
when(x=x.max)  Adoc1:x=0; 
when(x=x.min)  Inoi1:x=0;
when(x=x.max)  Inoi1:x=0; 
when(x=x.min)  Inoc1:x=0;
when(x=x.max)  Inoc1:x=0; 
// ADDITIONAL PARAMETERS
private real ZEROT = 0 mmol/cm^2;// Removes Initial Conditions from Parameter List
private real ZEROM = 0 mM;	// Removes Initial Conditions from Parameter List
//     ENZYME CONVERSION PARAMETERS
real Etot	= 0.001 mM,		// Enzyme Concentration in cell
     kf1	= 10 mM^(-1)*s^(-1),	// rate const ado+enz->ecmplx 
     kb1 	= 10 s^(-1),		// rate const ecmplx->ado+enz
     kf2 	= 10 sec^(-1),		// rate const ecmplx->ino+enz
     kb2 	= 0 mM^(-1)*s^(-1);	// rate const ino+enz->ecmplx
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
real Enz1(t,x) mM;		// Concentration of Enz1 in Vc
real ECmplx1(t,x) mM;		// Concentration of ECmplx1 in Vc
real TAdoi1(t,x) mmol/cm^2;		// Transporter complex
real TAdoc1(t,x) mmol/cm^2;		// Transporter complex
real TInoi1(t,x) mmol/cm^2;		// Transporter complex
real TInoc1(t,x) mmol/cm^2;		// Transporter complex
real Ti1(t,x) mmol/cm^2;		// Free transporter
real Tc1(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz1	= Etot;
when(t=t.min) ECmplx1	= ZEROM;
when(t=t.min) TAdoi1 	= ZEROT;
when(t=t.min) TAdoc1 	= ZEROT;
when(t=t.min) TInoi1 	= ZEROT;
when(t=t.min) TInoc1 	= ZEROT;
when(t=t.min) Ti1	= Ttot/2;
when(t=t.min) Tc1	= Ttot/2;
when(x=x.min)  Enz1:x=0;
when(x=x.max)  Enz1:x=0; 
when(x=x.min)  ECmplx1:x=0;
when(x=x.max)  ECmplx1:x=0; 
when(x=x.min)  TAdoi1:x=0;
when(x=x.max)  TAdoi1:x=0; 
when(x=x.min)  TAdoc1:x=0;
when(x=x.max)  TAdoc1:x=0; 
when(x=x.min)  TInoi1:x=0;
when(x=x.max)  TInoi1:x=0; 
when(x=x.min)  TInoc1:x=0;
when(x=x.max)  TInoc1:x=0; 
when(x=x.min)  Ti1:x=0;
when(x=x.max)  Ti1:x=0; 
when(x=x.min)  Tc1:x=0;
when(x=x.max)  Tc1:x=0; 
// PDE CALCULATIONS
Adop1:t = -(Fi1*L/Vp)*Adop1:x
        + DpAdo*Adop1:x:x
        +PSgAdo/Vp*(Adoi1-Adop1);
Inop1:t = -(Fi1*L/Vp)*Inop1:x
        + DpIno*Inop1:x:x
        +PSgIno/Vp*(Inoi1-Inop1);
Adoi1:t = DiAdo*Adoi1:x:x
       +PSgAdo/Vi*(Adop1-Adoi1)
       +PSpcAdo/Vc*(Adoc1-Adoi1)
       +(-konAdoi*Adoi1*Ti1 + kofAdoi*TAdoi1)*SoVi;
Adoc1:t = DcAdo*Adoc1:x:x
       +PSpcAdo/Vc*(Adoi1-Adoc1)
       -Gado2ino/Vc*Adoc1
       +(-konAdoc*Adoc1*Tc1 + kofAdoc*TAdoc1)*SoVc
       -kf1*Adoc1*Enz1 + kb1*ECmplx1;
Inoi1:t = DiIno*Inoi1:x:x
       +PSgIno/Vi*(Inop1-Inoi1)
       +PSpcIno/Vc*(Inoc1-Inoi1)
       +(-konInoi*Inoi1*Ti1 + kofInoi*TInoi1)*SoVi;
Inoc1:t = DcIno*Inoc1:x:x
       +PSpcIno/Vc*(Inoi1-Inoc1)
       +Gado2ino/Vc*Adoc1
       +(-konInoc*Inoc1*Tc1 + kofInoc*TInoc1)*SoVc
       +kf2*ECmplx1 - kb2*Inoc1*Enz1;
Ti1:t = -konAdoi*Adoi1*Ti1 + kofAdoi*TAdoi1 
     -konInoi*Inoi1*Ti1 + kofInoi*TInoi1 
     - kTi2c*Ti1 + kTc2i*Tc1;
TAdoi1:t = konAdoi*Adoi1*Ti1 - kofAdoi*TAdoi1
        - kTAdoi2c*TAdoi1 + kTAdoc2i*TAdoc1;
Tc1:t = -konAdoc*Adoc1*Tc1 + kofAdoc*TAdoc1 
     -konInoc*Inoc1*Tc1 + kofInoc*TInoc1 
     +kTi2c*Ti1 - kTc2i*Tc1;
TAdoc1:t = konAdoc*Adoc1*Tc1 - kofAdoc*TAdoc1
        +kTAdoi2c*TAdoi1 - kTAdoc2i*TAdoc1;
TInoi1:t = konInoi*Inoi1*Ti1 - kofInoi*TInoi1
        - kTInoi2c*TInoi1 + kTInoc2i*TInoc1;
TInoc1:t = konInoc*Inoc1*Tc1 - kofInoc*TInoc1
        +kTInoi2c*TInoi1 - kTInoc2i*TInoc1;
Enz1:t = -(kf1*Adoc1 + kb2*Inoc1)*Enz1
      + (kb1 + kf2)*ECmplx1;
ECmplx1:t = (kf1*Adoc1 + kb2*Inoc1)*Enz1
         - (kb1 + kf2)*ECmplx1;
real Adoout2(t) mM;		// Outflowing concentration 
real Inoout2(t) mM;		// Outflowing concentration 
real Adop2(t,x) mM;		// Concentration 
real Adoi2(t,x) mM;		// Concentration 
real Adoc2(t,x) mM;		// Concentration 
real Inop2(t,x) mM;		// Concentration 
real Inoi2(t,x) mM;		// Concentration 
real Inoc2(t,x) mM;		// Concentration 
when(t=t.min)  Adop2=0;
when(t=t.min)  Adoi2=0;
when(t=t.min)  Adoc2=0;
when(t=t.min)  Inop2=0;
when(t=t.min)  Inoi2=0;
when(t=t.min)  Inoc2=0;
when  (x=x.min) (-Fi2*L/Vp)*(Adop2-Adoin)+DpAdo*Adop2:x = 0; 
when  (x=x.max) { Adop2:x = 0; Adoout2 = Adop2;}
when  (x=x.min) (-Fi2*L/Vp)*(Inop2-Inoin)+DpIno*Inop2:x = 0; 
when  (x=x.max) { Inop2:x = 0; Inoout2 = Inop2;}
when(x=x.min)  Adoi2:x=0;
when(x=x.max)  Adoi2:x=0; 
when(x=x.min)  Adoc2:x=0;
when(x=x.max)  Adoc2:x=0; 
when(x=x.min)  Inoi2:x=0;
when(x=x.max)  Inoi2:x=0; 
when(x=x.min)  Inoc2:x=0;
when(x=x.max)  Inoc2:x=0; 
real Enz2(t,x) mM;		// Concentration of Enz2 in Vc
real ECmplx2(t,x) mM;		// Concentration of ECmplx2 in Vc
real TAdoi2(t,x) mmol/cm^2;		// Transporter complex
real TAdoc2(t,x) mmol/cm^2;		// Transporter complex
real TInoi2(t,x) mmol/cm^2;		// Transporter complex
real TInoc2(t,x) mmol/cm^2;		// Transporter complex
real Ti2(t,x) mmol/cm^2;		// Free transporter
real Tc2(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz2	= Etot;
when(t=t.min) ECmplx2	= ZEROM;
when(t=t.min) TAdoi2 	= ZEROT;
when(t=t.min) TAdoc2 	= ZEROT;
when(t=t.min) TInoi2 	= ZEROT;
when(t=t.min) TInoc2 	= ZEROT;
when(t=t.min) Ti2	= Ttot/2;
when(t=t.min) Tc2	= Ttot/2;
when(x=x.min)  Enz2:x=0;
when(x=x.max)  Enz2:x=0; 
when(x=x.min)  ECmplx2:x=0;
when(x=x.max)  ECmplx2:x=0; 
when(x=x.min)  TAdoi2:x=0;
when(x=x.max)  TAdoi2:x=0; 
when(x=x.min)  TAdoc2:x=0;
when(x=x.max)  TAdoc2:x=0; 
when(x=x.min)  TInoi2:x=0;
when(x=x.max)  TInoi2:x=0; 
when(x=x.min)  TInoc2:x=0;
when(x=x.max)  TInoc2:x=0; 
when(x=x.min)  Ti2:x=0;
when(x=x.max)  Ti2:x=0; 
when(x=x.min)  Tc2:x=0;
when(x=x.max)  Tc2:x=0; 
Adop2:t = -(Fi2*L/Vp)*Adop2:x
        + DpAdo*Adop2:x:x
        +PSgAdo/Vp*(Adoi2-Adop2);
Inop2:t = -(Fi2*L/Vp)*Inop2:x
        + DpIno*Inop2:x:x
        +PSgIno/Vp*(Inoi2-Inop2);
Adoi2:t = DiAdo*Adoi2:x:x
       +PSgAdo/Vi*(Adop2-Adoi2)
       +PSpcAdo/Vc*(Adoc2-Adoi2)
       +(-konAdoi*Adoi2*Ti2 + kofAdoi*TAdoi2)*SoVi;
Adoc2:t = DcAdo*Adoc2:x:x
       +PSpcAdo/Vc*(Adoi2-Adoc2)
       -Gado2ino/Vc*Adoc2
       +(-konAdoc*Adoc2*Tc2 + kofAdoc*TAdoc2)*SoVc
       -kf1*Adoc2*Enz2 + kb1*ECmplx2;
Inoi2:t = DiIno*Inoi2:x:x
       +PSgIno/Vi*(Inop2-Inoi2)
       +PSpcIno/Vc*(Inoc2-Inoi2)
       +(-konInoi*Inoi2*Ti2 + kofInoi*TInoi2)*SoVi;
Inoc2:t = DcIno*Inoc2:x:x
       +PSpcIno/Vc*(Inoi2-Inoc2)
       +Gado2ino/Vc*Adoc2
       +(-konInoc*Inoc2*Tc2 + kofInoc*TInoc2)*SoVc
       +kf2*ECmplx2 - kb2*Inoc2*Enz2;
Ti2:t = -konAdoi*Adoi2*Ti2 + kofAdoi*TAdoi2 
     -konInoi*Inoi2*Ti2 + kofInoi*TInoi2 
     - kTi2c*Ti2 + kTc2i*Tc2;
TAdoi2:t = konAdoi*Adoi2*Ti2 - kofAdoi*TAdoi2
        - kTAdoi2c*TAdoi2 + kTAdoc2i*TAdoc2;
Tc2:t = -konAdoc*Adoc2*Tc2 + kofAdoc*TAdoc2 
     -konInoc*Inoc2*Tc2 + kofInoc*TInoc2 
     +kTi2c*Ti2 - kTc2i*Tc2;
TAdoc2:t = konAdoc*Adoc2*Tc2 - kofAdoc*TAdoc2
        +kTAdoi2c*TAdoi2 - kTAdoc2i*TAdoc2;
TInoi2:t = konInoi*Inoi2*Ti2 - kofInoi*TInoi2
        - kTInoi2c*TInoi2 + kTInoc2i*TInoc2;
TInoc2:t = konInoc*Inoc2*Tc2 - kofInoc*TInoc2
        +kTInoi2c*TInoi2 - kTInoc2i*TInoc2;
Enz2:t = -(kf1*Adoc2 + kb2*Inoc2)*Enz2
      + (kb1 + kf2)*ECmplx2;
ECmplx2:t = (kf1*Adoc2 + kb2*Inoc2)*Enz2
         - (kb1 + kf2)*ECmplx2;
real Adoout3(t) mM;		// Outflowing concentration 
real Inoout3(t) mM;		// Outflowing concentration 
real Adop3(t,x) mM;		// Concentration 
real Adoi3(t,x) mM;		// Concentration 
real Adoc3(t,x) mM;		// Concentration 
real Inop3(t,x) mM;		// Concentration 
real Inoi3(t,x) mM;		// Concentration 
real Inoc3(t,x) mM;		// Concentration 
when(t=t.min)  Adop3=0;
when(t=t.min)  Adoi3=0;
when(t=t.min)  Adoc3=0;
when(t=t.min)  Inop3=0;
when(t=t.min)  Inoi3=0;
when(t=t.min)  Inoc3=0;
when  (x=x.min) (-Fi3*L/Vp)*(Adop3-Adoin)+DpAdo*Adop3:x = 0; 
when  (x=x.max) { Adop3:x = 0; Adoout3 = Adop3;}
when  (x=x.min) (-Fi3*L/Vp)*(Inop3-Inoin)+DpIno*Inop3:x = 0; 
when  (x=x.max) { Inop3:x = 0; Inoout3 = Inop3;}
when(x=x.min)  Adoi3:x=0;
when(x=x.max)  Adoi3:x=0; 
when(x=x.min)  Adoc3:x=0;
when(x=x.max)  Adoc3:x=0; 
when(x=x.min)  Inoi3:x=0;
when(x=x.max)  Inoi3:x=0; 
when(x=x.min)  Inoc3:x=0;
when(x=x.max)  Inoc3:x=0; 
real Enz3(t,x) mM;		// Concentration of Enz3 in Vc
real ECmplx3(t,x) mM;		// Concentration of ECmplx3 in Vc
real TAdoi3(t,x) mmol/cm^2;		// Transporter complex
real TAdoc3(t,x) mmol/cm^2;		// Transporter complex
real TInoi3(t,x) mmol/cm^2;		// Transporter complex
real TInoc3(t,x) mmol/cm^2;		// Transporter complex
real Ti3(t,x) mmol/cm^2;		// Free transporter
real Tc3(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz3	= Etot;
when(t=t.min) ECmplx3	= ZEROM;
when(t=t.min) TAdoi3 	= ZEROT;
when(t=t.min) TAdoc3 	= ZEROT;
when(t=t.min) TInoi3 	= ZEROT;
when(t=t.min) TInoc3 	= ZEROT;
when(t=t.min) Ti3	= Ttot/2;
when(t=t.min) Tc3	= Ttot/2;
when(x=x.min)  Enz3:x=0;
when(x=x.max)  Enz3:x=0; 
when(x=x.min)  ECmplx3:x=0;
when(x=x.max)  ECmplx3:x=0; 
when(x=x.min)  TAdoi3:x=0;
when(x=x.max)  TAdoi3:x=0; 
when(x=x.min)  TAdoc3:x=0;
when(x=x.max)  TAdoc3:x=0; 
when(x=x.min)  TInoi3:x=0;
when(x=x.max)  TInoi3:x=0; 
when(x=x.min)  TInoc3:x=0;
when(x=x.max)  TInoc3:x=0; 
when(x=x.min)  Ti3:x=0;
when(x=x.max)  Ti3:x=0; 
when(x=x.min)  Tc3:x=0;
when(x=x.max)  Tc3:x=0; 
Adop3:t = -(Fi3*L/Vp)*Adop3:x
        + DpAdo*Adop3:x:x
        +PSgAdo/Vp*(Adoi3-Adop3);
Inop3:t = -(Fi3*L/Vp)*Inop3:x
        + DpIno*Inop3:x:x
        +PSgIno/Vp*(Inoi3-Inop3);
Adoi3:t = DiAdo*Adoi3:x:x
       +PSgAdo/Vi*(Adop3-Adoi3)
       +PSpcAdo/Vc*(Adoc3-Adoi3)
       +(-konAdoi*Adoi3*Ti3 + kofAdoi*TAdoi3)*SoVi;
Adoc3:t = DcAdo*Adoc3:x:x
       +PSpcAdo/Vc*(Adoi3-Adoc3)
       -Gado2ino/Vc*Adoc3
       +(-konAdoc*Adoc3*Tc3 + kofAdoc*TAdoc3)*SoVc
       -kf1*Adoc3*Enz3 + kb1*ECmplx3;
Inoi3:t = DiIno*Inoi3:x:x
       +PSgIno/Vi*(Inop3-Inoi3)
       +PSpcIno/Vc*(Inoc3-Inoi3)
       +(-konInoi*Inoi3*Ti3 + kofInoi*TInoi3)*SoVi;
Inoc3:t = DcIno*Inoc3:x:x
       +PSpcIno/Vc*(Inoi3-Inoc3)
       +Gado2ino/Vc*Adoc3
       +(-konInoc*Inoc3*Tc3 + kofInoc*TInoc3)*SoVc
       +kf2*ECmplx3 - kb2*Inoc3*Enz3;
Ti3:t = -konAdoi*Adoi3*Ti3 + kofAdoi*TAdoi3 
     -konInoi*Inoi3*Ti3 + kofInoi*TInoi3 
     - kTi2c*Ti3 + kTc2i*Tc3;
TAdoi3:t = konAdoi*Adoi3*Ti3 - kofAdoi*TAdoi3
        - kTAdoi2c*TAdoi3 + kTAdoc2i*TAdoc3;
Tc3:t = -konAdoc*Adoc3*Tc3 + kofAdoc*TAdoc3 
     -konInoc*Inoc3*Tc3 + kofInoc*TInoc3 
     +kTi2c*Ti3 - kTc2i*Tc3;
TAdoc3:t = konAdoc*Adoc3*Tc3 - kofAdoc*TAdoc3
        +kTAdoi2c*TAdoi3 - kTAdoc2i*TAdoc3;
TInoi3:t = konInoi*Inoi3*Ti3 - kofInoi*TInoi3
        - kTInoi2c*TInoi3 + kTInoc2i*TInoc3;
TInoc3:t = konInoc*Inoc3*Tc3 - kofInoc*TInoc3
        +kTInoi2c*TInoi3 - kTInoc2i*TInoc3;
Enz3:t = -(kf1*Adoc3 + kb2*Inoc3)*Enz3
      + (kb1 + kf2)*ECmplx3;
ECmplx3:t = (kf1*Adoc3 + kb2*Inoc3)*Enz3
         - (kb1 + kf2)*ECmplx3;
real Adoout4(t) mM;		// Outflowing concentration 
real Inoout4(t) mM;		// Outflowing concentration 
real Adop4(t,x) mM;		// Concentration 
real Adoi4(t,x) mM;		// Concentration 
real Adoc4(t,x) mM;		// Concentration 
real Inop4(t,x) mM;		// Concentration 
real Inoi4(t,x) mM;		// Concentration 
real Inoc4(t,x) mM;		// Concentration 
when(t=t.min)  Adop4=0;
when(t=t.min)  Adoi4=0;
when(t=t.min)  Adoc4=0;
when(t=t.min)  Inop4=0;
when(t=t.min)  Inoi4=0;
when(t=t.min)  Inoc4=0;
when  (x=x.min) (-Fi4*L/Vp)*(Adop4-Adoin)+DpAdo*Adop4:x = 0; 
when  (x=x.max) { Adop4:x = 0; Adoout4 = Adop4;}
when  (x=x.min) (-Fi4*L/Vp)*(Inop4-Inoin)+DpIno*Inop4:x = 0; 
when  (x=x.max) { Inop4:x = 0; Inoout4 = Inop4;}
when(x=x.min)  Adoi4:x=0;
when(x=x.max)  Adoi4:x=0; 
when(x=x.min)  Adoc4:x=0;
when(x=x.max)  Adoc4:x=0; 
when(x=x.min)  Inoi4:x=0;
when(x=x.max)  Inoi4:x=0; 
when(x=x.min)  Inoc4:x=0;
when(x=x.max)  Inoc4:x=0; 
real Enz4(t,x) mM;		// Concentration of Enz4 in Vc
real ECmplx4(t,x) mM;		// Concentration of ECmplx4 in Vc
real TAdoi4(t,x) mmol/cm^2;		// Transporter complex
real TAdoc4(t,x) mmol/cm^2;		// Transporter complex
real TInoi4(t,x) mmol/cm^2;		// Transporter complex
real TInoc4(t,x) mmol/cm^2;		// Transporter complex
real Ti4(t,x) mmol/cm^2;		// Free transporter
real Tc4(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz4	= Etot;
when(t=t.min) ECmplx4	= ZEROM;
when(t=t.min) TAdoi4 	= ZEROT;
when(t=t.min) TAdoc4 	= ZEROT;
when(t=t.min) TInoi4 	= ZEROT;
when(t=t.min) TInoc4 	= ZEROT;
when(t=t.min) Ti4	= Ttot/2;
when(t=t.min) Tc4	= Ttot/2;
when(x=x.min)  Enz4:x=0;
when(x=x.max)  Enz4:x=0; 
when(x=x.min)  ECmplx4:x=0;
when(x=x.max)  ECmplx4:x=0; 
when(x=x.min)  TAdoi4:x=0;
when(x=x.max)  TAdoi4:x=0; 
when(x=x.min)  TAdoc4:x=0;
when(x=x.max)  TAdoc4:x=0; 
when(x=x.min)  TInoi4:x=0;
when(x=x.max)  TInoi4:x=0; 
when(x=x.min)  TInoc4:x=0;
when(x=x.max)  TInoc4:x=0; 
when(x=x.min)  Ti4:x=0;
when(x=x.max)  Ti4:x=0; 
when(x=x.min)  Tc4:x=0;
when(x=x.max)  Tc4:x=0; 
Adop4:t = -(Fi4*L/Vp)*Adop4:x
        + DpAdo*Adop4:x:x
        +PSgAdo/Vp*(Adoi4-Adop4);
Inop4:t = -(Fi4*L/Vp)*Inop4:x
        + DpIno*Inop4:x:x
        +PSgIno/Vp*(Inoi4-Inop4);
Adoi4:t = DiAdo*Adoi4:x:x
       +PSgAdo/Vi*(Adop4-Adoi4)
       +PSpcAdo/Vc*(Adoc4-Adoi4)
       +(-konAdoi*Adoi4*Ti4 + kofAdoi*TAdoi4)*SoVi;
Adoc4:t = DcAdo*Adoc4:x:x
       +PSpcAdo/Vc*(Adoi4-Adoc4)
       -Gado2ino/Vc*Adoc4
       +(-konAdoc*Adoc4*Tc4 + kofAdoc*TAdoc4)*SoVc
       -kf1*Adoc4*Enz4 + kb1*ECmplx4;
Inoi4:t = DiIno*Inoi4:x:x
       +PSgIno/Vi*(Inop4-Inoi4)
       +PSpcIno/Vc*(Inoc4-Inoi4)
       +(-konInoi*Inoi4*Ti4 + kofInoi*TInoi4)*SoVi;
Inoc4:t = DcIno*Inoc4:x:x
       +PSpcIno/Vc*(Inoi4-Inoc4)
       +Gado2ino/Vc*Adoc4
       +(-konInoc*Inoc4*Tc4 + kofInoc*TInoc4)*SoVc
       +kf2*ECmplx4 - kb2*Inoc4*Enz4;
Ti4:t = -konAdoi*Adoi4*Ti4 + kofAdoi*TAdoi4 
     -konInoi*Inoi4*Ti4 + kofInoi*TInoi4 
     - kTi2c*Ti4 + kTc2i*Tc4;
TAdoi4:t = konAdoi*Adoi4*Ti4 - kofAdoi*TAdoi4
        - kTAdoi2c*TAdoi4 + kTAdoc2i*TAdoc4;
Tc4:t = -konAdoc*Adoc4*Tc4 + kofAdoc*TAdoc4 
     -konInoc*Inoc4*Tc4 + kofInoc*TInoc4 
     +kTi2c*Ti4 - kTc2i*Tc4;
TAdoc4:t = konAdoc*Adoc4*Tc4 - kofAdoc*TAdoc4
        +kTAdoi2c*TAdoi4 - kTAdoc2i*TAdoc4;
TInoi4:t = konInoi*Inoi4*Ti4 - kofInoi*TInoi4
        - kTInoi2c*TInoi4 + kTInoc2i*TInoc4;
TInoc4:t = konInoc*Inoc4*Tc4 - kofInoc*TInoc4
        +kTInoi2c*TInoi4 - kTInoc2i*TInoc4;
Enz4:t = -(kf1*Adoc4 + kb2*Inoc4)*Enz4
      + (kb1 + kf2)*ECmplx4;
ECmplx4:t = (kf1*Adoc4 + kb2*Inoc4)*Enz4
         - (kb1 + kf2)*ECmplx4;
real Adoout5(t) mM;		// Outflowing concentration 
real Inoout5(t) mM;		// Outflowing concentration 
real Adop5(t,x) mM;		// Concentration 
real Adoi5(t,x) mM;		// Concentration 
real Adoc5(t,x) mM;		// Concentration 
real Inop5(t,x) mM;		// Concentration 
real Inoi5(t,x) mM;		// Concentration 
real Inoc5(t,x) mM;		// Concentration 
when(t=t.min)  Adop5=0;
when(t=t.min)  Adoi5=0;
when(t=t.min)  Adoc5=0;
when(t=t.min)  Inop5=0;
when(t=t.min)  Inoi5=0;
when(t=t.min)  Inoc5=0;
when  (x=x.min) (-Fi5*L/Vp)*(Adop5-Adoin)+DpAdo*Adop5:x = 0; 
when  (x=x.max) { Adop5:x = 0; Adoout5 = Adop5;}
when  (x=x.min) (-Fi5*L/Vp)*(Inop5-Inoin)+DpIno*Inop5:x = 0; 
when  (x=x.max) { Inop5:x = 0; Inoout5 = Inop5;}
when(x=x.min)  Adoi5:x=0;
when(x=x.max)  Adoi5:x=0; 
when(x=x.min)  Adoc5:x=0;
when(x=x.max)  Adoc5:x=0; 
when(x=x.min)  Inoi5:x=0;
when(x=x.max)  Inoi5:x=0; 
when(x=x.min)  Inoc5:x=0;
when(x=x.max)  Inoc5:x=0; 
real Enz5(t,x) mM;		// Concentration of Enz5 in Vc
real ECmplx5(t,x) mM;		// Concentration of ECmplx5 in Vc
real TAdoi5(t,x) mmol/cm^2;		// Transporter complex
real TAdoc5(t,x) mmol/cm^2;		// Transporter complex
real TInoi5(t,x) mmol/cm^2;		// Transporter complex
real TInoc5(t,x) mmol/cm^2;		// Transporter complex
real Ti5(t,x) mmol/cm^2;		// Free transporter
real Tc5(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz5	= Etot;
when(t=t.min) ECmplx5	= ZEROM;
when(t=t.min) TAdoi5 	= ZEROT;
when(t=t.min) TAdoc5 	= ZEROT;
when(t=t.min) TInoi5 	= ZEROT;
when(t=t.min) TInoc5 	= ZEROT;
when(t=t.min) Ti5	= Ttot/2;
when(t=t.min) Tc5	= Ttot/2;
when(x=x.min)  Enz5:x=0;
when(x=x.max)  Enz5:x=0; 
when(x=x.min)  ECmplx5:x=0;
when(x=x.max)  ECmplx5:x=0; 
when(x=x.min)  TAdoi5:x=0;
when(x=x.max)  TAdoi5:x=0; 
when(x=x.min)  TAdoc5:x=0;
when(x=x.max)  TAdoc5:x=0; 
when(x=x.min)  TInoi5:x=0;
when(x=x.max)  TInoi5:x=0; 
when(x=x.min)  TInoc5:x=0;
when(x=x.max)  TInoc5:x=0; 
when(x=x.min)  Ti5:x=0;
when(x=x.max)  Ti5:x=0; 
when(x=x.min)  Tc5:x=0;
when(x=x.max)  Tc5:x=0; 
Adop5:t = -(Fi5*L/Vp)*Adop5:x
        + DpAdo*Adop5:x:x
        +PSgAdo/Vp*(Adoi5-Adop5);
Inop5:t = -(Fi5*L/Vp)*Inop5:x
        + DpIno*Inop5:x:x
        +PSgIno/Vp*(Inoi5-Inop5);
Adoi5:t = DiAdo*Adoi5:x:x
       +PSgAdo/Vi*(Adop5-Adoi5)
       +PSpcAdo/Vc*(Adoc5-Adoi5)
       +(-konAdoi*Adoi5*Ti5 + kofAdoi*TAdoi5)*SoVi;
Adoc5:t = DcAdo*Adoc5:x:x
       +PSpcAdo/Vc*(Adoi5-Adoc5)
       -Gado2ino/Vc*Adoc5
       +(-konAdoc*Adoc5*Tc5 + kofAdoc*TAdoc5)*SoVc
       -kf1*Adoc5*Enz5 + kb1*ECmplx5;
Inoi5:t = DiIno*Inoi5:x:x
       +PSgIno/Vi*(Inop5-Inoi5)
       +PSpcIno/Vc*(Inoc5-Inoi5)
       +(-konInoi*Inoi5*Ti5 + kofInoi*TInoi5)*SoVi;
Inoc5:t = DcIno*Inoc5:x:x
       +PSpcIno/Vc*(Inoi5-Inoc5)
       +Gado2ino/Vc*Adoc5
       +(-konInoc*Inoc5*Tc5 + kofInoc*TInoc5)*SoVc
       +kf2*ECmplx5 - kb2*Inoc5*Enz5;
Ti5:t = -konAdoi*Adoi5*Ti5 + kofAdoi*TAdoi5 
     -konInoi*Inoi5*Ti5 + kofInoi*TInoi5 
     - kTi2c*Ti5 + kTc2i*Tc5;
TAdoi5:t = konAdoi*Adoi5*Ti5 - kofAdoi*TAdoi5
        - kTAdoi2c*TAdoi5 + kTAdoc2i*TAdoc5;
Tc5:t = -konAdoc*Adoc5*Tc5 + kofAdoc*TAdoc5 
     -konInoc*Inoc5*Tc5 + kofInoc*TInoc5 
     +kTi2c*Ti5 - kTc2i*Tc5;
TAdoc5:t = konAdoc*Adoc5*Tc5 - kofAdoc*TAdoc5
        +kTAdoi2c*TAdoi5 - kTAdoc2i*TAdoc5;
TInoi5:t = konInoi*Inoi5*Ti5 - kofInoi*TInoi5
        - kTInoi2c*TInoi5 + kTInoc2i*TInoc5;
TInoc5:t = konInoc*Inoc5*Tc5 - kofInoc*TInoc5
        +kTInoi2c*TInoi5 - kTInoc2i*TInoc5;
Enz5:t = -(kf1*Adoc5 + kb2*Inoc5)*Enz5
      + (kb1 + kf2)*ECmplx5;
ECmplx5:t = (kf1*Adoc5 + kb2*Inoc5)*Enz5
         - (kb1 + kf2)*ECmplx5;
real Adoout6(t) mM;		// Outflowing concentration 
real Inoout6(t) mM;		// Outflowing concentration 
real Adop6(t,x) mM;		// Concentration 
real Adoi6(t,x) mM;		// Concentration 
real Adoc6(t,x) mM;		// Concentration 
real Inop6(t,x) mM;		// Concentration 
real Inoi6(t,x) mM;		// Concentration 
real Inoc6(t,x) mM;		// Concentration 
when(t=t.min)  Adop6=0;
when(t=t.min)  Adoi6=0;
when(t=t.min)  Adoc6=0;
when(t=t.min)  Inop6=0;
when(t=t.min)  Inoi6=0;
when(t=t.min)  Inoc6=0;
when  (x=x.min) (-Fi6*L/Vp)*(Adop6-Adoin)+DpAdo*Adop6:x = 0; 
when  (x=x.max) { Adop6:x = 0; Adoout6 = Adop6;}
when  (x=x.min) (-Fi6*L/Vp)*(Inop6-Inoin)+DpIno*Inop6:x = 0; 
when  (x=x.max) { Inop6:x = 0; Inoout6 = Inop6;}
when(x=x.min)  Adoi6:x=0;
when(x=x.max)  Adoi6:x=0; 
when(x=x.min)  Adoc6:x=0;
when(x=x.max)  Adoc6:x=0; 
when(x=x.min)  Inoi6:x=0;
when(x=x.max)  Inoi6:x=0; 
when(x=x.min)  Inoc6:x=0;
when(x=x.max)  Inoc6:x=0; 
real Enz6(t,x) mM;		// Concentration of Enz6 in Vc
real ECmplx6(t,x) mM;		// Concentration of ECmplx6 in Vc
real TAdoi6(t,x) mmol/cm^2;		// Transporter complex
real TAdoc6(t,x) mmol/cm^2;		// Transporter complex
real TInoi6(t,x) mmol/cm^2;		// Transporter complex
real TInoc6(t,x) mmol/cm^2;		// Transporter complex
real Ti6(t,x) mmol/cm^2;		// Free transporter
real Tc6(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz6	= Etot;
when(t=t.min) ECmplx6	= ZEROM;
when(t=t.min) TAdoi6 	= ZEROT;
when(t=t.min) TAdoc6 	= ZEROT;
when(t=t.min) TInoi6 	= ZEROT;
when(t=t.min) TInoc6 	= ZEROT;
when(t=t.min) Ti6	= Ttot/2;
when(t=t.min) Tc6	= Ttot/2;
when(x=x.min)  Enz6:x=0;
when(x=x.max)  Enz6:x=0; 
when(x=x.min)  ECmplx6:x=0;
when(x=x.max)  ECmplx6:x=0; 
when(x=x.min)  TAdoi6:x=0;
when(x=x.max)  TAdoi6:x=0; 
when(x=x.min)  TAdoc6:x=0;
when(x=x.max)  TAdoc6:x=0; 
when(x=x.min)  TInoi6:x=0;
when(x=x.max)  TInoi6:x=0; 
when(x=x.min)  TInoc6:x=0;
when(x=x.max)  TInoc6:x=0; 
when(x=x.min)  Ti6:x=0;
when(x=x.max)  Ti6:x=0; 
when(x=x.min)  Tc6:x=0;
when(x=x.max)  Tc6:x=0; 
Adop6:t = -(Fi6*L/Vp)*Adop6:x
        + DpAdo*Adop6:x:x
        +PSgAdo/Vp*(Adoi6-Adop6);
Inop6:t = -(Fi6*L/Vp)*Inop6:x
        + DpIno*Inop6:x:x
        +PSgIno/Vp*(Inoi6-Inop6);
Adoi6:t = DiAdo*Adoi6:x:x
       +PSgAdo/Vi*(Adop6-Adoi6)
       +PSpcAdo/Vc*(Adoc6-Adoi6)
       +(-konAdoi*Adoi6*Ti6 + kofAdoi*TAdoi6)*SoVi;
Adoc6:t = DcAdo*Adoc6:x:x
       +PSpcAdo/Vc*(Adoi6-Adoc6)
       -Gado2ino/Vc*Adoc6
       +(-konAdoc*Adoc6*Tc6 + kofAdoc*TAdoc6)*SoVc
       -kf1*Adoc6*Enz6 + kb1*ECmplx6;
Inoi6:t = DiIno*Inoi6:x:x
       +PSgIno/Vi*(Inop6-Inoi6)
       +PSpcIno/Vc*(Inoc6-Inoi6)
       +(-konInoi*Inoi6*Ti6 + kofInoi*TInoi6)*SoVi;
Inoc6:t = DcIno*Inoc6:x:x
       +PSpcIno/Vc*(Inoi6-Inoc6)
       +Gado2ino/Vc*Adoc6
       +(-konInoc*Inoc6*Tc6 + kofInoc*TInoc6)*SoVc
       +kf2*ECmplx6 - kb2*Inoc6*Enz6;
Ti6:t = -konAdoi*Adoi6*Ti6 + kofAdoi*TAdoi6 
     -konInoi*Inoi6*Ti6 + kofInoi*TInoi6 
     - kTi2c*Ti6 + kTc2i*Tc6;
TAdoi6:t = konAdoi*Adoi6*Ti6 - kofAdoi*TAdoi6
        - kTAdoi2c*TAdoi6 + kTAdoc2i*TAdoc6;
Tc6:t = -konAdoc*Adoc6*Tc6 + kofAdoc*TAdoc6 
     -konInoc*Inoc6*Tc6 + kofInoc*TInoc6 
     +kTi2c*Ti6 - kTc2i*Tc6;
TAdoc6:t = konAdoc*Adoc6*Tc6 - kofAdoc*TAdoc6
        +kTAdoi2c*TAdoi6 - kTAdoc2i*TAdoc6;
TInoi6:t = konInoi*Inoi6*Ti6 - kofInoi*TInoi6
        - kTInoi2c*TInoi6 + kTInoc2i*TInoc6;
TInoc6:t = konInoc*Inoc6*Tc6 - kofInoc*TInoc6
        +kTInoi2c*TInoi6 - kTInoc2i*TInoc6;
Enz6:t = -(kf1*Adoc6 + kb2*Inoc6)*Enz6
      + (kb1 + kf2)*ECmplx6;
ECmplx6:t = (kf1*Adoc6 + kb2*Inoc6)*Enz6
         - (kb1 + kf2)*ECmplx6;
real Adoout7(t) mM;		// Outflowing concentration 
real Inoout7(t) mM;		// Outflowing concentration 
real Adop7(t,x) mM;		// Concentration 
real Adoi7(t,x) mM;		// Concentration 
real Adoc7(t,x) mM;		// Concentration 
real Inop7(t,x) mM;		// Concentration 
real Inoi7(t,x) mM;		// Concentration 
real Inoc7(t,x) mM;		// Concentration 
when(t=t.min)  Adop7=0;
when(t=t.min)  Adoi7=0;
when(t=t.min)  Adoc7=0;
when(t=t.min)  Inop7=0;
when(t=t.min)  Inoi7=0;
when(t=t.min)  Inoc7=0;
when  (x=x.min) (-Fi7*L/Vp)*(Adop7-Adoin)+DpAdo*Adop7:x = 0; 
when  (x=x.max) { Adop7:x = 0; Adoout7 = Adop7;}
when  (x=x.min) (-Fi7*L/Vp)*(Inop7-Inoin)+DpIno*Inop7:x = 0; 
when  (x=x.max) { Inop7:x = 0; Inoout7 = Inop7;}
when(x=x.min)  Adoi7:x=0;
when(x=x.max)  Adoi7:x=0; 
when(x=x.min)  Adoc7:x=0;
when(x=x.max)  Adoc7:x=0; 
when(x=x.min)  Inoi7:x=0;
when(x=x.max)  Inoi7:x=0; 
when(x=x.min)  Inoc7:x=0;
when(x=x.max)  Inoc7:x=0; 
real Enz7(t,x) mM;		// Concentration of Enz7 in Vc
real ECmplx7(t,x) mM;		// Concentration of ECmplx7 in Vc
real TAdoi7(t,x) mmol/cm^2;		// Transporter complex
real TAdoc7(t,x) mmol/cm^2;		// Transporter complex
real TInoi7(t,x) mmol/cm^2;		// Transporter complex
real TInoc7(t,x) mmol/cm^2;		// Transporter complex
real Ti7(t,x) mmol/cm^2;		// Free transporter
real Tc7(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz7	= Etot;
when(t=t.min) ECmplx7	= ZEROM;
when(t=t.min) TAdoi7 	= ZEROT;
when(t=t.min) TAdoc7 	= ZEROT;
when(t=t.min) TInoi7 	= ZEROT;
when(t=t.min) TInoc7 	= ZEROT;
when(t=t.min) Ti7	= Ttot/2;
when(t=t.min) Tc7	= Ttot/2;
when(x=x.min)  Enz7:x=0;
when(x=x.max)  Enz7:x=0; 
when(x=x.min)  ECmplx7:x=0;
when(x=x.max)  ECmplx7:x=0; 
when(x=x.min)  TAdoi7:x=0;
when(x=x.max)  TAdoi7:x=0; 
when(x=x.min)  TAdoc7:x=0;
when(x=x.max)  TAdoc7:x=0; 
when(x=x.min)  TInoi7:x=0;
when(x=x.max)  TInoi7:x=0; 
when(x=x.min)  TInoc7:x=0;
when(x=x.max)  TInoc7:x=0; 
when(x=x.min)  Ti7:x=0;
when(x=x.max)  Ti7:x=0; 
when(x=x.min)  Tc7:x=0;
when(x=x.max)  Tc7:x=0; 
Adop7:t = -(Fi7*L/Vp)*Adop7:x
        + DpAdo*Adop7:x:x
        +PSgAdo/Vp*(Adoi7-Adop7);
Inop7:t = -(Fi7*L/Vp)*Inop7:x
        + DpIno*Inop7:x:x
        +PSgIno/Vp*(Inoi7-Inop7);
Adoi7:t = DiAdo*Adoi7:x:x
       +PSgAdo/Vi*(Adop7-Adoi7)
       +PSpcAdo/Vc*(Adoc7-Adoi7)
       +(-konAdoi*Adoi7*Ti7 + kofAdoi*TAdoi7)*SoVi;
Adoc7:t = DcAdo*Adoc7:x:x
       +PSpcAdo/Vc*(Adoi7-Adoc7)
       -Gado2ino/Vc*Adoc7
       +(-konAdoc*Adoc7*Tc7 + kofAdoc*TAdoc7)*SoVc
       -kf1*Adoc7*Enz7 + kb1*ECmplx7;
Inoi7:t = DiIno*Inoi7:x:x
       +PSgIno/Vi*(Inop7-Inoi7)
       +PSpcIno/Vc*(Inoc7-Inoi7)
       +(-konInoi*Inoi7*Ti7 + kofInoi*TInoi7)*SoVi;
Inoc7:t = DcIno*Inoc7:x:x
       +PSpcIno/Vc*(Inoi7-Inoc7)
       +Gado2ino/Vc*Adoc7
       +(-konInoc*Inoc7*Tc7 + kofInoc*TInoc7)*SoVc
       +kf2*ECmplx7 - kb2*Inoc7*Enz7;
Ti7:t = -konAdoi*Adoi7*Ti7 + kofAdoi*TAdoi7 
     -konInoi*Inoi7*Ti7 + kofInoi*TInoi7 
     - kTi2c*Ti7 + kTc2i*Tc7;
TAdoi7:t = konAdoi*Adoi7*Ti7 - kofAdoi*TAdoi7
        - kTAdoi2c*TAdoi7 + kTAdoc2i*TAdoc7;
Tc7:t = -konAdoc*Adoc7*Tc7 + kofAdoc*TAdoc7 
     -konInoc*Inoc7*Tc7 + kofInoc*TInoc7 
     +kTi2c*Ti7 - kTc2i*Tc7;
TAdoc7:t = konAdoc*Adoc7*Tc7 - kofAdoc*TAdoc7
        +kTAdoi2c*TAdoi7 - kTAdoc2i*TAdoc7;
TInoi7:t = konInoi*Inoi7*Ti7 - kofInoi*TInoi7
        - kTInoi2c*TInoi7 + kTInoc2i*TInoc7;
TInoc7:t = konInoc*Inoc7*Tc7 - kofInoc*TInoc7
        +kTInoi2c*TInoi7 - kTInoc2i*TInoc7;
Enz7:t = -(kf1*Adoc7 + kb2*Inoc7)*Enz7
      + (kb1 + kf2)*ECmplx7;
ECmplx7:t = (kf1*Adoc7 + kb2*Inoc7)*Enz7
         - (kb1 + kf2)*ECmplx7;
real Adoout8(t) mM;		// Outflowing concentration 
real Inoout8(t) mM;		// Outflowing concentration 
real Adop8(t,x) mM;		// Concentration 
real Adoi8(t,x) mM;		// Concentration 
real Adoc8(t,x) mM;		// Concentration 
real Inop8(t,x) mM;		// Concentration 
real Inoi8(t,x) mM;		// Concentration 
real Inoc8(t,x) mM;		// Concentration 
when(t=t.min)  Adop8=0;
when(t=t.min)  Adoi8=0;
when(t=t.min)  Adoc8=0;
when(t=t.min)  Inop8=0;
when(t=t.min)  Inoi8=0;
when(t=t.min)  Inoc8=0;
when  (x=x.min) (-Fi8*L/Vp)*(Adop8-Adoin)+DpAdo*Adop8:x = 0; 
when  (x=x.max) { Adop8:x = 0; Adoout8 = Adop8;}
when  (x=x.min) (-Fi8*L/Vp)*(Inop8-Inoin)+DpIno*Inop8:x = 0; 
when  (x=x.max) { Inop8:x = 0; Inoout8 = Inop8;}
when(x=x.min)  Adoi8:x=0;
when(x=x.max)  Adoi8:x=0; 
when(x=x.min)  Adoc8:x=0;
when(x=x.max)  Adoc8:x=0; 
when(x=x.min)  Inoi8:x=0;
when(x=x.max)  Inoi8:x=0; 
when(x=x.min)  Inoc8:x=0;
when(x=x.max)  Inoc8:x=0; 
real Enz8(t,x) mM;		// Concentration of Enz8 in Vc
real ECmplx8(t,x) mM;		// Concentration of ECmplx8 in Vc
real TAdoi8(t,x) mmol/cm^2;		// Transporter complex
real TAdoc8(t,x) mmol/cm^2;		// Transporter complex
real TInoi8(t,x) mmol/cm^2;		// Transporter complex
real TInoc8(t,x) mmol/cm^2;		// Transporter complex
real Ti8(t,x) mmol/cm^2;		// Free transporter
real Tc8(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz8	= Etot;
when(t=t.min) ECmplx8	= ZEROM;
when(t=t.min) TAdoi8 	= ZEROT;
when(t=t.min) TAdoc8 	= ZEROT;
when(t=t.min) TInoi8 	= ZEROT;
when(t=t.min) TInoc8 	= ZEROT;
when(t=t.min) Ti8	= Ttot/2;
when(t=t.min) Tc8	= Ttot/2;
when(x=x.min)  Enz8:x=0;
when(x=x.max)  Enz8:x=0; 
when(x=x.min)  ECmplx8:x=0;
when(x=x.max)  ECmplx8:x=0; 
when(x=x.min)  TAdoi8:x=0;
when(x=x.max)  TAdoi8:x=0; 
when(x=x.min)  TAdoc8:x=0;
when(x=x.max)  TAdoc8:x=0; 
when(x=x.min)  TInoi8:x=0;
when(x=x.max)  TInoi8:x=0; 
when(x=x.min)  TInoc8:x=0;
when(x=x.max)  TInoc8:x=0; 
when(x=x.min)  Ti8:x=0;
when(x=x.max)  Ti8:x=0; 
when(x=x.min)  Tc8:x=0;
when(x=x.max)  Tc8:x=0; 
Adop8:t = -(Fi8*L/Vp)*Adop8:x
        + DpAdo*Adop8:x:x
        +PSgAdo/Vp*(Adoi8-Adop8);
Inop8:t = -(Fi8*L/Vp)*Inop8:x
        + DpIno*Inop8:x:x
        +PSgIno/Vp*(Inoi8-Inop8);
Adoi8:t = DiAdo*Adoi8:x:x
       +PSgAdo/Vi*(Adop8-Adoi8)
       +PSpcAdo/Vc*(Adoc8-Adoi8)
       +(-konAdoi*Adoi8*Ti8 + kofAdoi*TAdoi8)*SoVi;
Adoc8:t = DcAdo*Adoc8:x:x
       +PSpcAdo/Vc*(Adoi8-Adoc8)
       -Gado2ino/Vc*Adoc8
       +(-konAdoc*Adoc8*Tc8 + kofAdoc*TAdoc8)*SoVc
       -kf1*Adoc8*Enz8 + kb1*ECmplx8;
Inoi8:t = DiIno*Inoi8:x:x
       +PSgIno/Vi*(Inop8-Inoi8)
       +PSpcIno/Vc*(Inoc8-Inoi8)
       +(-konInoi*Inoi8*Ti8 + kofInoi*TInoi8)*SoVi;
Inoc8:t = DcIno*Inoc8:x:x
       +PSpcIno/Vc*(Inoi8-Inoc8)
       +Gado2ino/Vc*Adoc8
       +(-konInoc*Inoc8*Tc8 + kofInoc*TInoc8)*SoVc
       +kf2*ECmplx8 - kb2*Inoc8*Enz8;
Ti8:t = -konAdoi*Adoi8*Ti8 + kofAdoi*TAdoi8 
     -konInoi*Inoi8*Ti8 + kofInoi*TInoi8 
     - kTi2c*Ti8 + kTc2i*Tc8;
TAdoi8:t = konAdoi*Adoi8*Ti8 - kofAdoi*TAdoi8
        - kTAdoi2c*TAdoi8 + kTAdoc2i*TAdoc8;
Tc8:t = -konAdoc*Adoc8*Tc8 + kofAdoc*TAdoc8 
     -konInoc*Inoc8*Tc8 + kofInoc*TInoc8 
     +kTi2c*Ti8 - kTc2i*Tc8;
TAdoc8:t = konAdoc*Adoc8*Tc8 - kofAdoc*TAdoc8
        +kTAdoi2c*TAdoi8 - kTAdoc2i*TAdoc8;
TInoi8:t = konInoi*Inoi8*Ti8 - kofInoi*TInoi8
        - kTInoi2c*TInoi8 + kTInoc2i*TInoc8;
TInoc8:t = konInoc*Inoc8*Tc8 - kofInoc*TInoc8
        +kTInoi2c*TInoi8 - kTInoc2i*TInoc8;
Enz8:t = -(kf1*Adoc8 + kb2*Inoc8)*Enz8
      + (kb1 + kf2)*ECmplx8;
ECmplx8:t = (kf1*Adoc8 + kb2*Inoc8)*Enz8
         - (kb1 + kf2)*ECmplx8;
real Adoout9(t) mM;		// Outflowing concentration 
real Inoout9(t) mM;		// Outflowing concentration 
real Adop9(t,x) mM;		// Concentration 
real Adoi9(t,x) mM;		// Concentration 
real Adoc9(t,x) mM;		// Concentration 
real Inop9(t,x) mM;		// Concentration 
real Inoi9(t,x) mM;		// Concentration 
real Inoc9(t,x) mM;		// Concentration 
when(t=t.min)  Adop9=0;
when(t=t.min)  Adoi9=0;
when(t=t.min)  Adoc9=0;
when(t=t.min)  Inop9=0;
when(t=t.min)  Inoi9=0;
when(t=t.min)  Inoc9=0;
when  (x=x.min) (-Fi9*L/Vp)*(Adop9-Adoin)+DpAdo*Adop9:x = 0; 
when  (x=x.max) { Adop9:x = 0; Adoout9 = Adop9;}
when  (x=x.min) (-Fi9*L/Vp)*(Inop9-Inoin)+DpIno*Inop9:x = 0; 
when  (x=x.max) { Inop9:x = 0; Inoout9 = Inop9;}
when(x=x.min)  Adoi9:x=0;
when(x=x.max)  Adoi9:x=0; 
when(x=x.min)  Adoc9:x=0;
when(x=x.max)  Adoc9:x=0; 
when(x=x.min)  Inoi9:x=0;
when(x=x.max)  Inoi9:x=0; 
when(x=x.min)  Inoc9:x=0;
when(x=x.max)  Inoc9:x=0; 
real Enz9(t,x) mM;		// Concentration of Enz9 in Vc
real ECmplx9(t,x) mM;		// Concentration of ECmplx9 in Vc
real TAdoi9(t,x) mmol/cm^2;		// Transporter complex
real TAdoc9(t,x) mmol/cm^2;		// Transporter complex
real TInoi9(t,x) mmol/cm^2;		// Transporter complex
real TInoc9(t,x) mmol/cm^2;		// Transporter complex
real Ti9(t,x) mmol/cm^2;		// Free transporter
real Tc9(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz9	= Etot;
when(t=t.min) ECmplx9	= ZEROM;
when(t=t.min) TAdoi9 	= ZEROT;
when(t=t.min) TAdoc9 	= ZEROT;
when(t=t.min) TInoi9 	= ZEROT;
when(t=t.min) TInoc9 	= ZEROT;
when(t=t.min) Ti9	= Ttot/2;
when(t=t.min) Tc9	= Ttot/2;
when(x=x.min)  Enz9:x=0;
when(x=x.max)  Enz9:x=0; 
when(x=x.min)  ECmplx9:x=0;
when(x=x.max)  ECmplx9:x=0; 
when(x=x.min)  TAdoi9:x=0;
when(x=x.max)  TAdoi9:x=0; 
when(x=x.min)  TAdoc9:x=0;
when(x=x.max)  TAdoc9:x=0; 
when(x=x.min)  TInoi9:x=0;
when(x=x.max)  TInoi9:x=0; 
when(x=x.min)  TInoc9:x=0;
when(x=x.max)  TInoc9:x=0; 
when(x=x.min)  Ti9:x=0;
when(x=x.max)  Ti9:x=0; 
when(x=x.min)  Tc9:x=0;
when(x=x.max)  Tc9:x=0; 
Adop9:t = -(Fi9*L/Vp)*Adop9:x
        + DpAdo*Adop9:x:x
        +PSgAdo/Vp*(Adoi9-Adop9);
Inop9:t = -(Fi9*L/Vp)*Inop9:x
        + DpIno*Inop9:x:x
        +PSgIno/Vp*(Inoi9-Inop9);
Adoi9:t = DiAdo*Adoi9:x:x
       +PSgAdo/Vi*(Adop9-Adoi9)
       +PSpcAdo/Vc*(Adoc9-Adoi9)
       +(-konAdoi*Adoi9*Ti9 + kofAdoi*TAdoi9)*SoVi;
Adoc9:t = DcAdo*Adoc9:x:x
       +PSpcAdo/Vc*(Adoi9-Adoc9)
       -Gado2ino/Vc*Adoc9
       +(-konAdoc*Adoc9*Tc9 + kofAdoc*TAdoc9)*SoVc
       -kf1*Adoc9*Enz9 + kb1*ECmplx9;
Inoi9:t = DiIno*Inoi9:x:x
       +PSgIno/Vi*(Inop9-Inoi9)
       +PSpcIno/Vc*(Inoc9-Inoi9)
       +(-konInoi*Inoi9*Ti9 + kofInoi*TInoi9)*SoVi;
Inoc9:t = DcIno*Inoc9:x:x
       +PSpcIno/Vc*(Inoi9-Inoc9)
       +Gado2ino/Vc*Adoc9
       +(-konInoc*Inoc9*Tc9 + kofInoc*TInoc9)*SoVc
       +kf2*ECmplx9 - kb2*Inoc9*Enz9;
Ti9:t = -konAdoi*Adoi9*Ti9 + kofAdoi*TAdoi9 
     -konInoi*Inoi9*Ti9 + kofInoi*TInoi9 
     - kTi2c*Ti9 + kTc2i*Tc9;
TAdoi9:t = konAdoi*Adoi9*Ti9 - kofAdoi*TAdoi9
        - kTAdoi2c*TAdoi9 + kTAdoc2i*TAdoc9;
Tc9:t = -konAdoc*Adoc9*Tc9 + kofAdoc*TAdoc9 
     -konInoc*Inoc9*Tc9 + kofInoc*TInoc9 
     +kTi2c*Ti9 - kTc2i*Tc9;
TAdoc9:t = konAdoc*Adoc9*Tc9 - kofAdoc*TAdoc9
        +kTAdoi2c*TAdoi9 - kTAdoc2i*TAdoc9;
TInoi9:t = konInoi*Inoi9*Ti9 - kofInoi*TInoi9
        - kTInoi2c*TInoi9 + kTInoc2i*TInoc9;
TInoc9:t = konInoc*Inoc9*Tc9 - kofInoc*TInoc9
        +kTInoi2c*TInoi9 - kTInoc2i*TInoc9;
Enz9:t = -(kf1*Adoc9 + kb2*Inoc9)*Enz9
      + (kb1 + kf2)*ECmplx9;
ECmplx9:t = (kf1*Adoc9 + kb2*Inoc9)*Enz9
         - (kb1 + kf2)*ECmplx9;
real Adoout10(t) mM;		// Outflowing concentration 
real Inoout10(t) mM;		// Outflowing concentration 
real Adop10(t,x) mM;		// Concentration 
real Adoi10(t,x) mM;		// Concentration 
real Adoc10(t,x) mM;		// Concentration 
real Inop10(t,x) mM;		// Concentration 
real Inoi10(t,x) mM;		// Concentration 
real Inoc10(t,x) mM;		// Concentration 
when(t=t.min)  Adop10=0;
when(t=t.min)  Adoi10=0;
when(t=t.min)  Adoc10=0;
when(t=t.min)  Inop10=0;
when(t=t.min)  Inoi10=0;
when(t=t.min)  Inoc10=0;
when  (x=x.min) (-Fi10*L/Vp)*(Adop10-Adoin)+DpAdo*Adop10:x = 0; 
when  (x=x.max) { Adop10:x = 0; Adoout10 = Adop10;}
when  (x=x.min) (-Fi10*L/Vp)*(Inop10-Inoin)+DpIno*Inop10:x = 0; 
when  (x=x.max) { Inop10:x = 0; Inoout10 = Inop10;}
when(x=x.min)  Adoi10:x=0;
when(x=x.max)  Adoi10:x=0; 
when(x=x.min)  Adoc10:x=0;
when(x=x.max)  Adoc10:x=0; 
when(x=x.min)  Inoi10:x=0;
when(x=x.max)  Inoi10:x=0; 
when(x=x.min)  Inoc10:x=0;
when(x=x.max)  Inoc10:x=0; 
real Enz10(t,x) mM;		// Concentration of Enz10 in Vc
real ECmplx10(t,x) mM;		// Concentration of ECmplx10 in Vc
real TAdoi10(t,x) mmol/cm^2;		// Transporter complex
real TAdoc10(t,x) mmol/cm^2;		// Transporter complex
real TInoi10(t,x) mmol/cm^2;		// Transporter complex
real TInoc10(t,x) mmol/cm^2;		// Transporter complex
real Ti10(t,x) mmol/cm^2;		// Free transporter
real Tc10(t,x) mmol/cm^2;		// Free transporter
when(t=t.min) Enz10	= Etot;
when(t=t.min) ECmplx10	= ZEROM;
when(t=t.min) TAdoi10 	= ZEROT;
when(t=t.min) TAdoc10 	= ZEROT;
when(t=t.min) TInoi10 	= ZEROT;
when(t=t.min) TInoc10 	= ZEROT;
when(t=t.min) Ti10	= Ttot/2;
when(t=t.min) Tc10	= Ttot/2;
when(x=x.min)  Enz10:x=0;
when(x=x.max)  Enz10:x=0; 
when(x=x.min)  ECmplx10:x=0;
when(x=x.max)  ECmplx10:x=0; 
when(x=x.min)  TAdoi10:x=0;
when(x=x.max)  TAdoi10:x=0; 
when(x=x.min)  TAdoc10:x=0;
when(x=x.max)  TAdoc10:x=0; 
when(x=x.min)  TInoi10:x=0;
when(x=x.max)  TInoi10:x=0; 
when(x=x.min)  TInoc10:x=0;
when(x=x.max)  TInoc10:x=0; 
when(x=x.min)  Ti10:x=0;
when(x=x.max)  Ti10:x=0; 
when(x=x.min)  Tc10:x=0;
when(x=x.max)  Tc10:x=0; 
Adop10:t = -(Fi10*L/Vp)*Adop10:x
         + DpAdo*Adop10:x:x
         +PSgAdo/Vp*(Adoi10-Adop10);
Inop10:t = -(Fi10*L/Vp)*Inop10:x
         + DpIno*Inop10:x:x
         +PSgIno/Vp*(Inoi10-Inop10);
Adoi10:t = DiAdo*Adoi10:x:x
       +PSgAdo/Vi*(Adop10-Adoi10)
       +PSpcAdo/Vc*(Adoc10-Adoi10)
       +(-konAdoi*Adoi10*Ti10 + kofAdoi*TAdoi10)*SoVi;
Adoc10:t = DcAdo*Adoc10:x:x
       +PSpcAdo/Vc*(Adoi10-Adoc10)
       -Gado2ino/Vc*Adoc10
       +(-konAdoc*Adoc10*Tc10 + kofAdoc*TAdoc10)*SoVc
       -kf1*Adoc10*Enz10 + kb1*ECmplx10;
Inoi10:t = DiIno*Inoi10:x:x
       +PSgIno/Vi*(Inop10-Inoi10)
       +PSpcIno/Vc*(Inoc10-Inoi10)
       +(-konInoi*Inoi10*Ti10 + kofInoi*TInoi10)*SoVi;
Inoc10:t = DcIno*Inoc10:x:x
       +PSpcIno/Vc*(Inoi10-Inoc10)
       +Gado2ino/Vc*Adoc10
       +(-konInoc*Inoc10*Tc10 + kofInoc*TInoc10)*SoVc
       +kf2*ECmplx10 - kb2*Inoc10*Enz10;
Ti10:t = -konAdoi*Adoi10*Ti10 + kofAdoi*TAdoi10 
     -konInoi*Inoi10*Ti10 + kofInoi*TInoi10 
     - kTi2c*Ti10 + kTc2i*Tc10;
TAdoi10:t = konAdoi*Adoi10*Ti10 - kofAdoi*TAdoi10
        - kTAdoi2c*TAdoi10 + kTAdoc2i*TAdoc10;
Tc10:t = -konAdoc*Adoc10*Tc10 + kofAdoc*TAdoc10 
     -konInoc*Inoc10*Tc10 + kofInoc*TInoc10 
     +kTi2c*Ti10 - kTc2i*Tc10;
TAdoc10:t = konAdoc*Adoc10*Tc10 - kofAdoc*TAdoc10
        +kTAdoi2c*TAdoi10 - kTAdoc2i*TAdoc10;
TInoi10:t = konInoi*Inoi10*Ti10 - kofInoi*TInoi10
        - kTInoi2c*TInoi10 + kTInoc2i*TInoc10;
TInoc10:t = konInoc*Inoc10*Tc10 - kofInoc*TInoc10
        +kTInoi2c*TInoi10 - kTInoc2i*TInoc10;
Enz10:t = -(kf1*Adoc10 + kb2*Inoc10)*Enz10
      + (kb1 + kf2)*ECmplx10;
ECmplx10:t = (kf1*Adoc10 + kb2*Inoc10)*Enz10
         - (kb1 + kf2)*ECmplx10;
//
// COLLECT OUTFLOWS SUMMED BY WEIGHTS
real AdooutTot(t) mM;
real InooutTot(t) mM;
real Adooutw1(t) mM;
real Adooutw2(t) mM;
real Adooutw3(t) mM;
real Adooutw4(t) mM;
real Adooutw5(t) mM;
real Adooutw6(t) mM;
real Adooutw7(t) mM;
real Adooutw8(t) mM;
real Adooutw9(t) mM;
real Adooutw10(t) mM;
real Inooutw1(t) mM;
real Inooutw2(t) mM;
real Inooutw3(t) mM;
real Inooutw4(t) mM;
real Inooutw5(t) mM;
real Inooutw6(t) mM;
real Inooutw7(t) mM;
real Inooutw8(t) mM;
real Inooutw9(t) mM;
real Inooutw10(t) mM;
Adooutw1 = Adoout1*wts(1)*normf(1);
Adooutw2 = Adoout2*wts(2)*normf(2);
Adooutw3 = Adoout3*wts(3)*normf(3);
Adooutw4 = Adoout4*wts(4)*normf(4);
Adooutw5 = Adoout5*wts(5)*normf(5);
Adooutw6 = Adoout6*wts(6)*normf(6);
Adooutw7 = Adoout7*wts(7)*normf(7);
Adooutw8 = Adoout8*wts(8)*normf(8);
Adooutw9 = Adoout9*wts(9)*normf(9);
Adooutw10 = Adoout10*wts(10)*normf(10);
Inooutw1 = Inoout1*wts(1)*normf(1);
Inooutw2 = Inoout2*wts(2)*normf(2);
Inooutw3 = Inoout3*wts(3)*normf(3);
Inooutw4 = Inoout4*wts(4)*normf(4);
Inooutw5 = Inoout5*wts(5)*normf(5);
Inooutw6 = Inoout6*wts(6)*normf(6);
Inooutw7 = Inoout7*wts(7)*normf(7);
Inooutw8 = Inoout8*wts(8)*normf(8);
Inooutw9 = Inoout9*wts(9)*normf(9);
Inooutw10 = Inoout10*wts(10)*normf(10);
AdooutTot = Adooutw1
          +Adooutw2
          +Adooutw3
          +Adooutw4
          +Adooutw5
          +Adooutw6
          +Adooutw7
          +Adooutw8
          +Adooutw9
          +Adooutw10;
InooutTot = Inooutw1
          +Inooutw2
          +Inooutw3
          +Inooutw4
          +Inooutw5
          +Inooutw6
          +Inooutw7
          +Inooutw8
          +Inooutw9
          +Inooutw10;
/*
   aiAdo  ( aoAdo)          area of input (output) curve
   tiAdo  ( to,  tsysAdo)   transit time of input (output, system)
   RDiAdo (RDoAdo, RDsysAdo)   relative dispersion of input (output, system) 
*/
real aiAdo, tiAdo,RDiAdo,aoAdo,toAdo,RDoAdo,tsysAdo,RDsysAdo;
curveStatistics(Adoin@t,AdooutTot@t,aiAdo,tiAdo,RDiAdo,aoAdo,toAdo,RDoAdo,tsysAdo,RDsysAdo);
/*
   aiIno  ( aoIno)          area of input (output) curve
   tiIno  ( to,  tsysIno)   transit time of input (output, system)
   RDiIno (RDoIno, RDsysIno)   relative dispersion of input (output, system) 
*/
real aiIno, tiIno,RDiIno,aoIno,toIno,RDoIno,tsysIno,RDsysIno;
curveStatistics(Inoin@t,InooutTot@t,aiIno,tiIno,RDiIno,aoIno,toIno,RDoIno,tsysIno,RDsysIno);
}
// This MML file generated from MultiflowAdo2Ino.mpc using MPC v1.01.
  
