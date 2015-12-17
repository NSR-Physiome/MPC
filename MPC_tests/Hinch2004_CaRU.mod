// MODEL NAME: Hinch2004_CuRU

// imported from SemGen

//%START units
import nsrunit;
unit conversion on;
unit mJ_per_mole_K=.001 kilogram^1*meter^2*second^(-2)*kelvin^(-1)*mole^(-1);
unit ms=.001 second^1;
unit um3_per_ms=1E-15 meter^3*second^(-1);
unit per_ms=1E3 second^(-1);
unit per_ms3=1E9 second^(-3);
unit mM=1 meter^(-3)*mole^1;
unit um3=1E-18 meter^3;
unit C_per_mole=1 second^1*ampere^1*mole^(-1);
unit um3_mM_per_ms=1E-15 second^(-1)*mole^1;
unit mM_per_ms=1E3 meter^(-3)*second^(-1)*mole^1;
//%END units

math main {

//%START timeDomain
	realDomain t ms; t.min=0; extern t.max; extern t.delta;
//%END timeDomain

//%START universalConsts
    real Cm = 1 uF/cm^2;  // membrane cap per area
    real F C_per_mole;    // Faraday constant
    F=96487;
    real F_Ca = 2*F;      // Faraday const, Ca2+ 
    real R mJ_per_mole_K; // Universal gas constant
    R=8314.5; 
//%END universalConsts

//%START externalInputs     
    real T kelvin;
    T=295;
    extern real Ca_o mM;   // Extracellular calcium (2+) concentration
    extern real Ca_i mM;   // Intracellular calcium (2+) concentration
    extern real Ca_SR mM;  // Calcium (2+) concentration in sarcoplasmic reticulum
    extern real minusI(t) uA/cm^2;
    extern real Vclamp(t) mV;            
    real clamp_0no_1yes = 0; 
    real Vdepolar = -90 mV;   
    real VV(t) mV, 
          V(t) mV;   // displacement of the cell membrane potential from rest (depolarization negative);
    when (t=t.min) { VV = if (clamp_0no_1yes=0) Vdepolar  else Vclamp; }
    V   = if (clamp_0no_1yes = 0) VV else Vclamp;
    real FVRT(t) dimensionless;    // Convenience variable that is equal to F*V/(R*T)
    real FVRT_Ca(t) dimensionless; // Convenience variable that is equal to (F*V/(R*T))*2
 
//%END externalInputs	

//%START CellGeometry
real Cap_mem = 1.534e-4 uF;  // Total membrane capacitance
real Area_mem um^2;          // membrane area
     Area_mem = Cap_mem/Cm;
real V_SR = 2.098e3 um^3;    // Vol of SR
real V_myo = 25.84e3 um^3;   // myocyte vol
      
//%END CellGeometry

//%START CaRU
	real J_L um3_per_ms;
	J_L=9.13e-4;
	real J_R um3_per_ms;
	J_R=0.02;
	real I_RyR(t) mM_per_ms;
	real g_D um3_per_ms;
	g_D=0.065;
	real I_LCC(t) mM_per_ms; // Total whole-cell current through the L-type calcium channels.
	real N dimensionless;
	N=50000;
	real expVL(t) dimensionless; // Convenience variable that is equal to exp((V-V_L)/del_VL)
	real c dimensionless; // Biasing to make inactivation a function of calcium 2+ in the dyadic space.
	c=0.01;
	real beta_poc(t) per_ms; //  Transition rate for three-state model of the ryanodine receptor 
                              //  when L-type calcium channel is open and ryanodine receptor is closed
	real t_L ms;
	t_L=1;
	real tau_R ms;
	tau_R=2.43;
	
	real alpha_p(t) per_ms;  // Transition rate for three-state model of the L-type calcium channel
	
	real mu_pcc per_ms;
	
	real a dimensionless; // Biasing to make inactivation a function of membrane voltage
	a=0.0625;
	
	real beta_m per_ms;   // Transition rate for three-state model of the ryanodine receptor
	real C_co(t) mM;  // Concentration of calcium (2+) in dyadic space when L-type calcium channel state is closed 
                       // and when ryanodine receptor state is open.
	real C_oc(t) mM;  // Concentration of calcium (2+) in the dyadic space when the L-type calcium channel state 
                       // is open and the ryanodine receptor state is closed
	real b dimensionless; // Biasing to make inactivation a function of membrane voltage
	b=14;
	real epsilon_pcc(t) per_ms;	
	real epsilon_pco(t) per_ms;
	real epsilon_m(t) per_ms;

	real V_L mV;       // Potential when half LCC open
	V_L=-2;	
	real del_VL mV;    // Width of opening potentials of L-type calcium channels
	del_VL=7;
	
	real phi_L dimensionless;
	phi_L=2.35;
	
	real mu_mcc per_ms;
	
	real phi_R dimensionless;
	phi_R=0.05;
	
	real mu_poc(t) per_ms;
	
	real theta_R dimensionless;
	theta_R=0.012;
	
	real tau_L ms;
	tau_L=650;
	
	real d dimensionless; // Biasing to make inactivation a function of calcium 2+ in the dyadic space
	d=100;
	
	real mu_moc(t) per_ms;
	real t_R ms;
	real alpha_m per_ms;   // Transition rate for three-state model of the L-type calcium channel
	real beta_pcc(t) per_ms;  // Transition rate for three-state model of the ryanodine receptor when 
                               // both the L-type calcium channel and ryanodine receptor are closed
	
	real K_L mM;
	K_L=0.22e-3;
	real K_RyR mM;
	K_RyR=41e-3;
	
	real J_Loc(t) um3_mM_per_ms;
	real y_oc(t) dimensionless;
	real J_L2(t) um3_mM_per_ms;
	real J_Loo(t) um3_mM_per_ms;
	real y_oo(t) dimensionless;
	
	real z_1(t) dimensionless;
	when(t=t.min) z_1=0.98859;
	
	real z_2(t) dimensionless;
	when(t=t.min) z_2=0.0087302;
	
	real J_L1(t) um3_mM_per_ms;
	real y_co(t) dimensionless;
	real y_cc(t) dimensionless;
	real r_3(t) per_ms;
	real r_7(t) per_ms;
	
	real z_3(t) dimensionless;
	when(t=t.min) z_3=0.0026566;
	
	real r_1(t) per_ms;
	real r_4 per_ms;
	real r_8(t) per_ms;
	real z_4(t) dimensionless;
	real r_6(t) per_ms;
	real r_5(t) per_ms;
	real r_2(t) per_ms;
	real denom(t) per_ms3;
	real J_Roo(t) um3_mM_per_ms;
	real J_Rco(t) um3_mM_per_ms;
	real J_R1(t) um3_mM_per_ms;
	real J_R3(t) um3_mM_per_ms;

	// <component name="membrane">
	FVRT_Ca=(2*FVRT);
	FVRT=(F*V/(R*T));

	// <component name="CaRU_Transitions">
	alpha_m=(phi_L/t_L);
	mu_mcc=(theta_R*d*(Ca_i^2+c*K_RyR^2)/(tau_R*(d*Ca_i^2+c*K_RyR^2)));
	expVL=exp((V-V_L)/del_VL);
	mu_pcc=((Ca_i^2+c*K_RyR^2)/(tau_R*(Ca_i^2+K_RyR^2)));
	t_R=(1.17*t_L);
	beta_pcc=(Ca_i^2/(t_R*(Ca_i^2+K_RyR^2)));
	mu_moc=(theta_R*d*(C_oc^2+c*K_RyR^2)/(tau_R*(d*C_oc^2+c*K_RyR^2)));
	beta_m=(phi_R/t_R);
	epsilon_pco=(C_co*(expVL+a)/(tau_L*K_L*(expVL+1)));
	mu_poc=((C_oc^2+c*K_RyR^2)/(tau_R*(C_oc^2+K_RyR^2)));
	beta_poc=(C_oc^2/(t_R*(C_oc^2+K_RyR^2)));
	alpha_p=(expVL/(t_L*(expVL+1)));
	epsilon_m=(b*(expVL+a)/(tau_L*(b*expVL+a)));
	epsilon_pcc=(Ca_i*(expVL+a)/(tau_L*K_L*(expVL+1)));

	// <component name="LCC_current">
	J_L2=(J_Loc*alpha_p/(alpha_p+alpha_m));
	I_LCC=((z_1*J_L1+z_2*J_L2)*N/V_myo);
	J_L1=(J_Loo*y_oo+J_Loc*y_oc);

	// <component name="DS_Calcium_Concentrations">
	C_oc=(if (abs(FVRT_Ca)>1E-9) (Ca_i+J_L/g_D*Ca_o*FVRT_Ca*exp((-1)*FVRT_Ca)/(1-exp((-1)*FVRT_Ca)))/(1+J_L/g_D*FVRT_Ca/(1-exp((-1)*FVRT_Ca))) else (Ca_i+J_L/g_D*Ca_o)/(1+J_L/g_D));
	C_co=((Ca_i+J_R/g_D*Ca_SR)/(1+J_R/g_D));

	// <component name="CaRU_reduced_states">
	z_3:t=(r_5*z_1-(r_6+r_3)*z_3+r_4*z_4);
	z_4=(1-z_1-z_2-z_3);
	r_3=(beta_m*mu_pcc/(beta_m+beta_pcc));
	z_1:t=((-1)*(r_1+r_5)*z_1+r_2*z_2+r_6*z_3);
	r_6=epsilon_m;
	r_5=(y_co*epsilon_pco+y_cc*epsilon_pcc);
	r_7=(alpha_m*epsilon_pcc/(alpha_p+alpha_m));
	r_4=mu_mcc;
	r_2=((alpha_p*mu_moc+alpha_m*mu_mcc)/(alpha_p+alpha_m));
	z_2:t=(r_1*z_1-(r_2+r_7)*z_2+r_8*z_4);
	r_8=epsilon_m;
	r_1=(y_oc*mu_poc+y_cc*mu_pcc);

	// <component name="CaRU_states">
	y_co=(alpha_m*(beta_pcc*(alpha_m+beta_m+beta_poc)+beta_poc*alpha_p)/denom);
	y_oc=(alpha_p*beta_m*(alpha_p+alpha_m+beta_m+beta_pcc)/denom);
	denom=((alpha_p+alpha_m)*((alpha_m+beta_m+beta_poc)*(beta_m+beta_pcc)+alpha_p*(beta_m+beta_poc)));
	y_oo=(alpha_p*(beta_poc*(alpha_p+beta_m+beta_pcc)+beta_pcc*alpha_m)/denom);
	y_cc=(alpha_m*beta_m*(alpha_m+alpha_p+beta_m+beta_poc)/denom);

	// <component name="LCC_and_RyR_fluxes">
	J_Loo=(if (abs(FVRT_Ca)>1E-5) J_L*FVRT_Ca/(1-exp((-1)*FVRT_Ca))*(Ca_o*exp((-1)*FVRT_Ca)-Ca_i+J_R/g_D*(Ca_o*exp((-1)*FVRT_Ca)-Ca_SR))/(1+J_R/g_D+J_L/g_D*FVRT_Ca/(1-exp(FVRT_Ca))) else J_L*1E-5/(1-exp((-1)*1E-5))*(Ca_o*exp((-1)*1E-5)-Ca_i+J_R/g_D*(Ca_o*exp((-1)*1E-5)-Ca_SR))/(1+J_R/g_D+J_L/g_D*1E-5/(1-exp((-1)*1E-5))));
	J_Loc=(if (abs(FVRT_Ca)>1E-5) J_L*FVRT_Ca/(1-exp((-1)*FVRT_Ca))*(Ca_o*exp((-1)*FVRT_Ca)-Ca_i)/(1+J_L/g_D*FVRT_Ca/(1-exp((-1)*FVRT_Ca))) else J_L*1E-5/(1-exp((-1)*1E-5))*(Ca_o*exp((-1)*1E-5)-Ca_i)/(1+J_L/g_D*1E-5/(1-exp((-1)*1E-5))));
	J_Roo=(if (abs(FVRT_Ca)>1E-5) J_R*(Ca_SR-Ca_i+J_L/g_D*FVRT_Ca/(1-exp((-1)*FVRT_Ca))*(Ca_SR-Ca_o*exp((-1)*FVRT_Ca)))/(1+J_R/g_D+J_L/g_D*FVRT_Ca/(1-exp((-1)*FVRT_Ca))) else J_R*(Ca_SR-Ca_i+J_L/g_D*1E-5/(1-exp((-1)*1E-5))*(Ca_SR-Ca_o*exp((-1)*1E-5)))/(1+J_R/g_D+J_L/g_D*1E-5/(1-exp((-1)*1E-5))));
	J_Rco=(J_R*(Ca_SR-Ca_i)/(1+J_R/g_D));

	// <component name="RyR_current">
	I_RyR=((z_1*J_R1+z_3*J_R3)*N/V_myo);
	J_R3=(J_Rco*beta_pcc/(beta_m+beta_pcc));
	J_R1=(y_oo*J_Roo+J_Rco*y_co);

//%END CaRU

// Needed for stand alone model:
  real Iion(t) uA/cm^2;
       Iion  = (I_RyR+I_LCC)*F_Ca*V_myo/Area_mem; 
  VV:t = if (clamp_0no_1yes = 0) (-minusI-Iion)/Cm else 0; // universal var

}

/*
----------------------------------------------
Codeword, Units, Value (if static), Definition
----------------------------------------------

a, dimensionless, 0.0625, Biasing to make inactivation a function of membrane voltage
alpha_m, per_ms, , Transition rate for three-state model of the L-type calcium channel
alpha_p, per_ms, , Transition rate for three-state model of the L-type calcium channel
b, dimensionless, 14, Biasing to make inactivation a function of membrane voltage
beta_m, per_ms, , Transition rate for three-state model of the ryanodine receptor
beta_pcc, per_ms, , Transition rate for three-state model of the ryanodine receptor when both the L-type calcium channel and ryanodine receptor are closed
beta_poc, per_ms, , Transition rate for three-state model of the ryanodine receptor when L-type calcium channel is open and ryanodine receptor is closed
c, dimensionless, 0.01, Biasing to make inactivation a function of calcium 2+ in the dyadic space
Ca_i, mM, , 
Ca_o, mM, , Extracellular calcium (2+) concentration
Ca_SR, mM, , 
C_co, mM, , Concentration of calcium (2+) in dyadic space when L-type calcium channel state is closed and when ryanodine receptor state is open
C_oc, mM, , Concentration of calcium (2+) in the dyadic space when the L-type calcium channel state is open and the ryanodine receptor state is closed
d, dimensionless, 100, Biasing to make inactivation a function of calcium 2+ in the dyadic space
del_VL, mV, 7, Width of opening potentials of L-type calcium channels
denom, per_ms3, , Convenience term for computing ycc, yco, yoc, and yoo variables
epsilon_m, per_ms, , Transition rate for three-state model of the L-type calcium channel
epsilon_pcc, per_ms, , Transition rate for three-state model of the L-type calcium channel
epsilon_pco, per_ms, , Transition rate for three-state model of the L-type calcium channel
expVL, dimensionless, , Convenience variable that is equal to exp((V-V_L)/del_VL)
F, C_per_mole, 96487, Faraday constant
FVRT, dimensionless, , Convenience variable that is equal to F*V/(R*T)
FVRT_Ca, dimensionless, , Convenience variable that is equal to (F*V/(R*T))*2. Also equivalent to [greek delta]*V in the source manuscript
g_D, um3_per_ms, 0.065, Calcium (2+) flux rate from dyadic space to cytosol
I_LCC, mM_per_ms, , 
I_RyR, mM_per_ms, , Total whole-cell current through the ryanodine receptors
J_L, um3_per_ms, 9.13e-4, Permeability of a single L-type calcium channel
J_L1, um3_mM_per_ms, , Convenience term for expressing equation for I_LCC?
J_L2, um3_mM_per_ms, , Convenience term for expressing equation for I_LCC?
J_Loc, um3_mM_per_ms, , Current through the L-type calcium channel when L-type calcium channel state is open and ryanodine receptor state is closed
J_Loo, um3_mM_per_ms, , Current through the L-type calcium channel when both the L-type calcium channel and ryanodine receptor are in the open state
J_R, um3_per_ms, 0.02, Permeability of a single ryanodine receptor
J_R1, um3_mM_per_ms, , Convenience term for expressing I_RyR equation?
J_R3, um3_mM_per_ms, , Convenience term for expressing I_RyR equation?
J_Rco, um3_mM_per_ms, , Current through the ryanodine receptor when L-type calcium channel state is closed and ryanodine receptor state is open
J_Roo, um3_mM_per_ms, , Current through the ryanodine receptor when both the L-type calcium channel and ryanodine receptor are in the open state
K_L, mM, 0.22e-3, Concentration of L-type calcium channels at inactivation
K_RyR, mM, 41e-3, Half concentration of activation for ryanodine receptors
mu_mcc, per_ms, , Transition rate for three-state model of the ryanodine receptor
mu_moc, per_ms, , Transition rate for three-state model of the ryanodine receptor
mu_pcc, per_ms, , Transition rate for three-state model of the ryanodine receptor
mu_poc, per_ms, , Transition rate for three-state model of the ryanodine receptor
N, dimensionless, , Number of calcium release units
phi_L, dimensionless, 2.35, Proportion of time closed in open mode for L-type calcium channels
phi_R, dimensionless, 0.05, Proportion of time closed in open mode for L-type calcium channels
R, mJ_per_mole_K, 8314.5, 
r_1, per_ms, , Rate constant for transition of calcium release unit from state 1 to state 2
r_2, per_ms, , Rate constant for transition of calcium release unit from state 2 to state 1
r_3, per_ms, , Rate constant for transition of calcium release unit from state 3 to state 4
r_4, per_ms, , Rate constant for transition of calcium release unit from state 4 to state 3
r_5, per_ms, , Rate constant for transition of calcium release unit from state 1 to state 3
r_6, per_ms, , Rate constant for transition of calcium release unit from state 3 to state 1
r_7, per_ms, , Rate constant for transition of calcium release unit from state 2 to state 4
r_8, per_ms, , Rate constant for transition of calcium release unit from state 4 to state 2
T, kelvin, 295, Ambient temperature
t_L, ms, 1, Time switching between closed and open states for L-type calcium channels
t_R, ms, , Time switching between closed and open states for ryanodine receptors
tau_L, ms, 650, Inactivation time for L-type calcium channels
tau_R, ms, 2.43, Inactivation time for ryanodine receptors
theta_R, dimensionless, 0.012, Reciprocal of proportion of time inactivated (closed) in open mode for ryanodine receptors
time, ms, , Simulation time domain
V, mV, , 
V_myo, um3, , 
V_L, mV, -2, 
y_cc, dimensionless, , The state of the calcium (2+) release unit when the L-type calcium channel and the ryanodine receptor are both in the closed state
y_co, dimensionless, , The state of the calcium (2+) release unit when the L-type calcium channel is in the closed state and the ryanodine receptor is in the open state
y_oc, dimensionless, , The state of the calcium (2+) release unit when the L-type calcium channel is in the open state and the ryanodine receptor is in the closed state
y_oo, dimensionless, , The state of the calcium (2+) release unit when both the L-type calcium channel and the ryanodine receptor are in the open state
z_1, dimensionless, , State 1 of calcium release unit: z1 = ycc union yco union yoc union yoo
z_2, dimensionless, , State 2 of calcium release channel: z2 = yoi union yci
z_3, dimensionless, , State 3 of calcium release unit: z3 = yic union yio
z_4, dimensionless, , State 4 of calcium release unit: z4 = yii
*/

