math ShortCodeLibrary {
real version = 1.0; 
/*--------------- BEGIN CONSTRUCT MODULAR PROGRAM CONSTRUCTOR LIBRARY
//------------------------------ HEADER
//%START           header
import nsrunit; unit conversion on;
math model {
//%END             header
//------------------------------ ODE DOMAINS
//%START           odeDomains
realDomain t s; t.min=0; t.max=30; t.delta = 0.1;
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
//------------------------------ REACTION A->B
//%START      reactionCalc
A:t = -G/V*A;
B:t = G/V*A;
//%END        reactionCalc
/*--------------- END CONSTRUCT MODULAR PROGRAM CONSTRUCTOR LIBRARY*/
} 
// This MML file generated from ShortCodeLibrary.mpc using MPC.
