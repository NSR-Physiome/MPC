//%START curveStatJava
//*****************************************************************************
/* Statistics on a curve */

source procedure curveStatistics(input@t,output@t; ai,ti,RDi,ao,tao,RDo,tsys,RDsys) {

/* ai  (ao)       area of input (output) curve
   ti  (to)       transit time of input (output) curve
   RDi (RDo)      relative dispersion of input (output) curve
   tsys           transit time of the system
   RDsys          relative dispersion of the system

   NOTA BENE, the results are not dimensionally scaled or checked.
   A basic assumption is that both the input and output have the
   same units.  

       language="java";
       maincode={{
           RegularGridData t = (RegularGridData) input.grid(0);
           double S1=0; double S2=0; double S3=0;
           double T1=0; double T2=0; double T3=0;
           for (int i=0; i<t.ct(); i++) {
               double time = t.min()+((double)i)*t.delta();
               double time2 = time*time;
               double timeO=time;
               double timeO2=timeO*timeO;
               double in=input.realVal(i);
               double ot=output.realVal(i);
               S1+=in;
               S2+=in*time;
               S3+=in*time2;
               T1+=ot;
               T2+=ot*timeO;
               T3+=ot*timeO2;
           }
           ai.set(t.delta()*S1);
           double tI=S2/S1;
           ti.set(tI);
           double rdI=( Math.sqrt((Math.abs(S3/S2)/(S2/S1))-1));
           RDi.set(rdI);

           ao.set(t.delta()*T1);
           double tO=T2/T1;
           tao.set(tO);
           double rdO=( Math.sqrt((Math.abs(T3/T2)/(T2/T1))-1));
           RDo.set(rdO);

           tsys.set(tO-tI);
           RDsys.set( Math.sqrt(rdO*tO*rdO*tO-rdI*tI*rdI*tI)/(tO-tI) );
         }};
} // END
//%END curveStatJava
import nsrunit;
unit conversion on;
math transit {
realDomain t sec; t.min=0; t.max=100; t.delta=0.1;
real C(t) mM;
real V=0.05 ml; 
real F= 1 ml/min;
when(t=t.min) {C=0;}
extern real Cin(t) mM; 
C:t=F*(Cin-C)/V;
real Cout(t) mM;
Cout=C;
//%START transitCalc
/*
   ai  ( ao)          area of input (output) curve
   ti  ( to,  tsys)   transit time of input (output, system)
   RDi (RDo, RDsys)   relative dispersion of input (output, system) 
*/
real ai, ti,RDi,ao,tao,RDo,tsys,RDsys;
curveStatistics(Cin@t,Cout@t,ai,ti,RDi,ao,tao,RDo,tsys,RDsys);
//%END transitCalc
}
