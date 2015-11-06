/* SHORT DESCRIPTION: Use A3comp to make adenosine to
   inosine model for plasma, isf, cell.
*/
import nsrunit; unit conversion on;
math Ado2Ino3comp {
//%START Ado2Ino3compModel
realDomain t s; t.min=0; t.max=30; t.delta = 0.1;
// PARAMETERS
real Flow = 1 ml/(g*min);	// Flow rate 
real PSgAdo = 3 ml/(g*min);	// Exchg rate 
real PSpcAdo = 5 ml/(g*min);	// Exchg rate 
real Vp = 0.05 ml/g;		// Volume of Vp
real Vi = 0.15 ml/g;		// Volume of Vi
real Vc = 0.60 ml/g;		// Volume of Vc
extern real Adoin(t) mM;	// Inflowing concentration
// DEPENDENT VARIABLES
real Adoout(t) mM;		// Outflowing concentration 
real Adop(t) mM;		// Concentration 
real Adoi(t) mM;		// Concentration 
real Adoc(t) mM;		// Concentration 
// INITIAL CONDITIONS
when(t=t.min)  Adop=0;
when(t=t.min)  Adoi=0;
when(t=t.min)  Adoc=0;
real PSgIno = 3 ml/(g*min);	// Exchg rate 
real PSpcIno = 5 ml/(g*min);	// Exchg rate 
extern real Inoin(t) mM;	// Inflowing concentration
real Inoout(t) mM;		// Outflowing concentration 
real Inop(t) mM;		// Concentration 
real Inoi(t) mM;		// Concentration 
real Inoc(t) mM;		// Concentration 
when(t=t.min)  Inop=0;
when(t=t.min)  Inoi=0;
when(t=t.min)  Inoc=0;
real Gc = 15 ml/(g*min);
//-----------------------------------------------------------------------------
// ODE CALCULATIONS
Adoout =Adop;
Adop:t = (Flow/Vp)*(Adoin-Adop)
       +PSgAdo/Vp*(Adoi-Adop);
Adoi:t = PSgAdo/Vi*(Adop-Adoi)
       +PSpcAdo/Vi*(Adoc-Adoi);
Adoc:t = PSpcAdo/Vc*(Adoi-Adoc)
         -Gc/Vc*Adoc;
Inoout =Inop;
Inop:t = (Flow/Vp)*(Inoin-Inop)
       +PSgIno/Vp*(Inoi-Inop);
Inoi:t = PSgIno/Vi*(Inop-Inoi)
       +PSpcIno/Vi*(Inoc-Inoi);
Inoc:t = PSpcIno/Vc*(Inoi-Inoc)
         +Gc/Vc*Adoc;
//%END Ado2Ino3compModel 
}
// This MML file generated from Ado2Ino3comp.mpc using MPC.
