/* SHORT DESCRIPTION: A two species (A and B) model in
3 regions (1=flowing, 2 and 3 stagnant) with A->B
in region 3 */
import nsrunit; unit conversion on;
math A2B { 
// INDEPENDENT VARIABLES
realDomain  t s; t.min=0; t.max=30; t.delta = 0.1;
real L = 0.1 cm;
real Ndivx = 31;
realDomain x cm ; x.min=0; x.max=L; x.ct = Ndivx;
//%START a2bParmsVars
// PARAMETERS
real Flow = 1 ml/(g*min); // Flow rate 
real PSA12 = 3 ml/(g*min); // Exchg rate 
real PSA23 = 5 ml/(g*min); // Exchg rate 
real PSB12 = 3 ml/(g*min); // Exchg rate 
real PSB23 = 5 ml/(g*min); // Exchg rate 
real Ga2b = 10 ml/(g*min); // Conversion rate 
real V1 = 0.05 ml/g; // Volume of V1
real V2 = 0.15 ml/g; // Volume of V2
real V3 = 0.60 ml/g; // Volume of V3
real DA1 = 1e-6 cm^2/sec; // Diffusion coeff 
real DA2 = 1e-6 cm^2/sec; // Diffusion coeff 
real DA3 = 1e-6 cm^2/sec; // Diffusion coeff 
real DB1 = 1e-6 cm^2/sec; // Diffusion coeff 
real DB2 = 1e-6 cm^2/sec; // Diffusion coeff 
real DB3 = 1e-6 cm^2/sec; // Diffusion coeff 
extern real Ain(t) mM; // Inflowing concentration
extern real Bin(t) mM; // Inflowing concentration
// DEPENDENT VARIABLES
real Aout(t) mM; // Outflowing concentration 
real Bout(t) mM; // Outflowing concentration 
real A1(t,x) mM; // Concentration 
real A2(t,x) mM; // Concentration 
real A3(t,x) mM; // Concentration 
real B1(t,x) mM; // Concentration 
real B2(t,x) mM; // Concentration 
real B3(t,x) mM; // Concentration 
// INITIAL CONDITIONS
when(t=t.min) A1=0;
when(t=t.min) A2=0;
when(t=t.min) A3=0;
when(t=t.min) B1=0;
when(t=t.min) B2=0;
when(t=t.min) B3=0;
// BOUNDARY CONDITIONS
when  (x=x.min) (-Flow*L/V1)*(A1-Ain)+DA1*A1:x = 0; 
when  (x=x.max) { A1:x = 0; Aout = A1;}
when  (x=x.min) (-Flow*L/V1)*(B1-Bin)+DB1*B1:x = 0; 
when  (x=x.max) { B1:x = 0; Bout = B1;}
when(x=x.min)  A2:x=0;
when(x=x.max)  A2:x=0; 
when(x=x.min)  A3:x=0;
when(x=x.max)  A3:x=0; 
when(x=x.min)  B2:x=0;
when(x=x.max)  B2:x=0; 
when(x=x.min)  B3:x=0;
when(x=x.max)  B3:x=0; 
//%END a2bParmsVars
//%START a2bCalc
// PDE CALCULATIONS
A1:t = -(Flow*L/V1)*A1:x
     + DA1*A1:x:x 
       +PSA12/V1*(A2-A1);
A2:t = DA2*A2:x:x 
       +PSA12/V2*(A1-A2)
       +PSA23/V3*(A3-A2);
B1:t = -(Flow*L/V1)*B1:x
     + DB1*B1:x:x 
       +PSB12/V1*(B2-B1);
B2:t = DB2*B2:x:x 
       +PSB12/V2*(B1-B2)
       +PSB23/V3*(B3-B2);
A3:t = DA3*A3:x:x 
       +PSA23/V3*(A2-A3)
       -Ga2b/V3*A3;
B3:t = DB3*B3:x:x 
       +PSB23/V3*(B2-B3)
       +Ga2b/V3*A3;
//%END a2bCalc
}
// This MML file generated from A2B.mpc using MPC v1.01.
  
