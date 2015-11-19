math ShortCodeLibrary {
real version = 1.0; 
//--------------- BEGIN CONSTRUCT MODULAR PROGRAM CONSTRUCTOR LIBRARY
//------------------------------ HEADER
//%START           header
import nsrunit; unit conversion on;
math model {
//%END             header
//------------------------------ ODE DOMAINS
//%START           odeDomains
realDomain t s; t.min=0; t.max=16; t.delta = 0.1;
//%END             odeDomains
//------------------------------ flowCalc
//%START           flowCalc
C:t = (F/V)*(Cin-C);
//%END             flowCalc
//------------------------------ EXCHANGE CACULATIONS
//%START      exchangeCalc
C1:t	= PS/V1*(C2-C1);
C2:t	= PS/V2*(C1-C2);
//%END        exchangeCalc
//
//----------- PDE calcs ----------
//------------------------------ PDE DOMAINS
//%START pdeDomains 
realDomain  t s; t.min=0; t.max=30; t.delta = 0.1; 
real L = 0.1 cm; 
real Ndivx = 31; 
realDomain x cm ; x.min=0; x.max=L; x.ct = Ndivx; 
//%END pdeDomains 
// ----------------------------- PDE FLOW BCs
//%START flowBC 
when  (x=x.min) (-F*L/V)*(C-Cin)+D*C:x = 0; 
when  (x=x.max) { C:x = 0; Cout = C;} 
//%END flowBC 
//------------------------------ PDE FLOW DIFFUSION calcs 
//%START flowDiffCalc 
C:t = D*C:x:x 
      -(F*L/V)*C:x; 
//%END flowDiffCalc 
//------------------------------ REACTION A->B
//%START      reactionCalc
real G = 5 ml/(g*min);     // Const reaction rate.
A:t = -G/V*A;
B:t = G/V*A;
//%END        reactionCalc
//----------- PDE calcs END ---------

//------------------------------ MM REACTION A->B
//%START      MMreactionCalc
real KmA =1.0 mM, VmaxA =2 umol/(g*min); // MM constant and max velocity of rxn
real G(t) ml/(g*min);                  // MM reaction rate 
G = (VmaxA/(KmA+A));  
A:t = -G*(A)/V;
B:t =  G*(A)/V;
//%END        MMreactionCalc

/*--------------- END CONSTRUCT MODULAR PROGRAM CONSTRUCTOR LIBRARY*/
} 

