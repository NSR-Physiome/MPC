/* MODEL NUMBER:
   MODEL NAME: FlowHet
   SHORT DESCRIPTION: Generates relative flows and weights for heterogeneous
   models. Used by Modular Program Constructor (MPC).
*/

import nsrunit; unit conversion on;

math FlowHet{

real PATHS = 20; // SET PATHS in presource file

//%START multiFlowCalc1
//
/* NOTA BENE!************************************
   To use this code in a MultiFlow model, you
   must remove the real declaration for flow in
   the model (e.g. remove real Fp = 1 ml/(g*min);
   from a BTEX10 when making MultiFlowBTEX10.
*/
   private real NPATH   = PATHS*1;
   public real NPATHS   = NPATH;
   real relFmin = 0.05 ;
   real relFmax = 3.0 ;
//
   realDomain NP ; NP.min=1; NP.max=NPATH; NP.delta=1;
   private NP.min, NP.max, NP.delta, NP.ct;
   private real userF(NP), userWt(NP);
//%END multiFlowCalc1
// USER DEFINED FLOWS (from MPC)
//%START multiFlowCalc2
//
   realDomain NPp1 ; NPp1.min=1; NPp1.max=NPATH+1; NPp1.delta=1;
   private NPp1.min, NPp1.max, NPp1.delta, NPp1.ct;
//
   realDomain rflo ; rflo.min=relFmin; rflo.max=relFmax; rflo.delta=0.01;
   private rflo.min, rflo.max,  rflo.delta, rflo.ct;
//
   extern real relFpdf(rflo) ;	// pdf for relative flows
   private real cumpdf(rflo) ;
   when(rflo=rflo.min) cumpdf=0;
   cumpdf:rflo=relFpdf;
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
   private real ends(NPp1) ;         // ends of each flow bin
   ends(NPp1) = if(spacing=1) exp(srelFmin+(NPp1-1)*swidth) else
             if(spacing=2) 1/ (srelFmin+(NPp1-1)*swidth) else
                           (srelFmin+(NPp1-1)*swidth);
//
// MAKE EACH RELATIVE FLOW THE CENTER OF ITS FLOW BIN
   real f(NP) ;       // relative flows
   f(NP) = if(spacing<>4) (ends(NP+1)+ends(NP))/2 else userF(NP);
//
// INTERPOLATE cumpdf(rflo) TO GET cumValues(NPp1);
   private real cumValues(NPp1) ;    // cumpdf(rflo)->cumValues(PATHS+1)
   cumValues(NPp1) = cumpdf(ends(NPp1));
//
// GENERATE WT FOR EACH BIN so that integral wt*df=1
   private real wt(NP) ;             // d(cumValues)/d(ends)
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
//%END multiFlowCalc2
}
/* DETAILED DESCRIPTION
// 
// ALGORITHM FOR HETEROGENEOUS FLOWS
   1.) An algorithm for heterogeneous flow:
   A probability density function (PDF) is generated. Suggested
   choices are LagNormal or Gaussian, with area = 1, tMean =1,
   and RD = .4 to .55. The PDF represents the distribution of
   the input concentration among the various paths.
//
   2.) A range of relative flows is selected, relFmin, and relFmax.
   Reasonable values are relFmin=0.05; relFmax=3.0. The user
   chooses the number of flow paths, PATHS. Reasonable values
   are 5 to 20. The reader chooses how the relative flows are
   generated: (RECOMMENDED) equally spaced in the logarithm of the 
   relative flows; equally spaced in the transit time of the relative flows
   (emphasizes the slowest speeds);  or equally spaced in the relative flows
   (preserves the shape of the PDF but ignores the slowest flows). 
//
   3.) The OUTPUT for the user is wts(NPATH), and normf(NPATH).
//
   4.) This program is designed to function with the MODULAR PROGRAM
   CONSTRUCTOR. The .mpc file using this code must be modified and 
   run everytime the number of relative flow paths is changed. 
//
   5.) Individual path flows should be summed as
       Couti*wti*normfi
*/

