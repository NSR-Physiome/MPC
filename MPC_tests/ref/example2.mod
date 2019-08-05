
import nsrunit; unit conversion on; // Use SI units
math example {                      // model declaration
// INDEPENDENT VARIABLES 
realDomain t s; t.min=0; t.max=16; t.delta = 0.1;
//%START a2bParmsVars       // Specify parameters and variables section
// PARAMETERS 
real Flow = 1 ml/(g*min);	    // Flow rate 
real PSA12 = 6 ml/(g*min);   // Conductances: PSA12,PSB12,PSC12
real PSB12 = 5 ml/(g*min);   // Conductances: PSA12,PSB12,PSC12
real PSC12 = 4 ml/(g*min);   // Conductances: PSA12,PSB12,PSC12
real V1 = 0.05 ml/g;		    // Volume of V1, V2 
real V2 = 0.05 ml/g;		    // Volume of V1, V2 
extern real Ain(t) mM;	    // Inflowing concentrations 
extern real Bin(t) mM;	    // Inflowing concentrations 
extern real Cin(t) mM;	    // Inflowing concentrations 
// DEPENDENT VARIABLES 
real A1(t) mM;		    // A1,A2,B1,B2,C1,C2 
real A2(t) mM;		    // A1,A2,B1,B2,C1,C2 
real B1(t) mM;		    // A1,A2,B1,B2,C1,C2 
real B2(t) mM;		    // A1,A2,B1,B2,C1,C2 
real C1(t) mM;		    // A1,A2,B1,B2,C1,C2 
real C2(t) mM;		    // A1,A2,B1,B2,C1,C2 
// INITIAL CONDITIONS (IC's) 
when(t=t.min)  A1=0;         // Defines IC's for the ODEs
when(t=t.min)  A2=0;         // Defines IC's for the ODEs
when(t=t.min)  B1=0;         // Defines IC's for the ODEs
when(t=t.min)  B2=0;         // Defines IC's for the ODEs
when(t=t.min)  C1=0;         // Defines IC's for the ODEs
when(t=t.min)  C2=0;         // Defines IC's for the ODEs
//%END a2bParmsVars         // End parameters and variables section
//%START a2bCalc            // Specify calculations section
// ODE CALCULATIONS 
A1:t = (Flow/V1)*(Ain-A1)
       +PSA12/V1*(A2-A1);
B1:t = (Flow/V1)*(Bin-B1)
       +PSB12/V1*(B2-B1);
C1:t = (Flow/V1)*(Cin-C1)
       +PSC12/V1*(C2-C1);
real Ga2b = 5 ml/(g*min);     // Const reaction rate.
A2:t = PSA12/V2*(A1-A2)
       -Ga2b/V2*A2;
real KmB2 =1.0 mM, VmaxB2 =2 umol/(g*min); // MM constant and max velocity of rxn
real Gb2c(t) ml/(g*min);                  // MM reaction rate 
Gb2c = (VmaxB2/(KmB2+B2));  
B2:t = PSB12/V2*(B1-B2)
       +Ga2b/V2*A2
       -Gb2c*(B2)/V2;
C2:t = PSC12/V2*(C1-C2)
       +Gb2c*(B2)/V2;
//%END a2bCalc
}  // curly bracket ends model

// This MML file generated from example2.mpc using MPC v1.02.
