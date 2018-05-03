/*
   MODEL NAME: FlowHet
   SHORT DESCRIPTION: Generates relative flows and weights for heterogeneous
   models. Used by Modular Program Constructor (MPC).
*/
//%START main
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
       "NUMBER=7" for 7 flows
   7.) To get the weighted flow in each path
ChemOutWeightedPathNumber = ChemOutPathNumber*wts(PathNumber)*normf(PathNumber)
ChemOutTot = ChemOutWeightedPathnumber

MPC's COLLECT(ChemoutTot)  will collect the total outflow

*/
import nsrunit; unit conversion on;
math Multiflow{

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
//%END main
}
// This MML file generated from Multiflow.mpc using MPC v1.01.
  
