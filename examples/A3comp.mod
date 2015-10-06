// SHORT DESCRIPTION: A single species ODE model in 3 
// compartments (1=flowing, 2 and 3 stagnant).
import nsrunit; unit conversion on;
math A3comp {
// INDEPENDENT VARIABLES
realDomain t s; t.min=0; t.max=30; t.delta = 0.1;
//%START A3compParmsVars
// PARAMETERS
real Flow = 1 ml/(g*min);	// Flow rate 
real PSA12 = 3 ml/(g*min);	// Exchg rate 
real PSA23 = 5 ml/(g*min);	// Exchg rate 
real Vol1 = 0.05 ml/g;		// Volume of Vol1
real Vol2 = 0.15 ml/g;		// Volume of Vol2
real Vol3 = 0.60 ml/g;		// Volume of Vol3
extern real Ain(t) mM;	// Inflowing concentration
// DEPENDENT VARIABLES
real Aout(t) mM;		// Outflowing concentration 
real A1(t) mM;		// Concentration 
real A2(t) mM;		// Concentration 
real A3(t) mM;		// Concentration 
// INITIAL CONDITIONS
when(t=t.min)  A1=0;
when(t=t.min)  A2=0;
when(t=t.min)  A3=0;
//%END A3compParmsVars
//%START A3compCalc
// ODE CALCULATIONS
Aout =A1;
A1:t = (Flow/Vol1)*(Ain-A1)
       +PSA12/Vol1*(A2-A1);
A2:t = PSA12/Vol2*(A1-A2)
       +PSA23/Vol2*(A3-A2);
A3:t = PSA23/Vol3*(A2-A3);
//%END A3compCalc
}
// This MML file generated from A3comp.mpc using MPC.
  
