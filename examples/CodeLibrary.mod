math ShortCodeLibrary {
real version = 1.0; 
/*--------------- BEGIN CONSTRUCT MODULAR PROGRAM CONSTRUCTOR LIBRARY
//------------------------------ PDE DOMAINS
//%START           pdeDomains
// INDEPENDENT VARIABLES
realDomain  t s; t.min=0; t.max=30; t.delta = 0.1;
real L = 0.1 cm;
real Ndivx = 31;
realDomain x cm ; x.min=0; x.max=L; x.ct = Ndivx;
//%END             pdeDomains
//------------------------------ BOUNDARY CONDITIONS
//%START           flowBC
when  (x=x.min) (-F*L/V)*(C-Cin)+D*C:x = 0; 
when  (x=x.max) { C:x = 0; Cout = C;}
//%END             flowBC

//%START       noFlowBC
when(x=x.min)  C:x=0;
when(x=x.max)  C:x=0; 
//%END         noFlowBC
//------------------------------ FLOW DIFF CALCULATION
//%START           flowDiffCalc
C:t =  -(F*L/V)*C:x
       + D*C:x:x ;
//%END             flowDiffCalc
//------------------------------ DIFFUSION CALCULATION
//%START       diffusionCalc
C:t = D*C:x:x ;
//%END         diffusionCalc
//------------------------------ EXCHANGE CACULATIONS
//%START      exchangeCalc
C1:t	= PS/V1*(C2-C1);
C2:t	= PS/V2*(C1-C2);
//%END        exchangeCalc
//------------------------------ CONSUME CALCULATION
//%START      consumeCalc
C:t = -(G/V)*C;
//%END        consumeCalc
//------------------------------ REACTION A->B
//%START      reactionCalc
A:t = -G/V*A;
B:t = G/V*A;
//%END        reactionCalc
//------------------------------ ON OFF MEMBRANE BINDING SITE
//%START           onOffMembraneCalc
M:t = (-kon*M*B + kof*MB)*SoV;
B:t =  -kon*M*B + kof*MB ;
MB:t =  kon*M*B - kof*MB;
//%END             onOffMembraneCalc
//------------------------------ CONFORMATIONAL CHANGE (FLIP)
//%START                  flipa2bCalc
a:t = - ka2b*a + kb2a*b;
b:t =   ka2b*a - kb2a*b;
//%END                    flipa2bCalc
//------------------------------ ENZYME CONVERSION
//%START    enzymeParms
real Etot = 0.1 mM;
real kf1 = 2 mM^(-1)*s^(-1);
real kb1 = 1 s^(-1);
real kf2 = 1 s^(-1);
real kb2 = 1 mM^(-1)*s^(-1);
//%END      enzymeParms
//%START    enzymeCalc
A:t = -kf1*A*Enzyme + kb1*Complex;
B:t =  kf2*Complex - kb2*B*Enzyme;
Enzyme:t = -(kf1*A + kb2*B)*Enzyme 
           + (kb1 + kf2)*Complex;
Complex:t = (kf1*A + kb2*B)*Enzyme 
           - (kb1 + kf2)*Complex;
//%END      enzymeCalc
/*--------------- END CONSTRUCT MODULAR PROGRAM CONSTRUCTOR LIBRARY
} 
// This MML file generated from ShortCodeLibrary.mpc using MPC.
