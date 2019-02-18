%---------------------------------------------------------------------------------------------------%
%---                                  Explicit Euler Method                                     ----%
%---------------------------------------------------------------------------------------------------%
function [V]=Vexata(a,b,c)
global  m n h G_Na G_K G_L Ena Ek El nt dt V C_M Iext

%-----                                Approximate solution of V                                -----%
d=dt/C_M;
for i=1:nt-1
AlphaMv(i)=( (25-V(i))/10 )/( exp( (25-V(i))/10 )-1 );  BetaMv(i) =4*exp(-V(i)/18);
AlphaNv(i)=(10-V(i))/( 100*( exp( (10-V(i))/10 )-1) );  BetaNv(i) =0.125*exp(-V(i)/80);
AlphaHv(i)=0.07*exp( -V(i)/20 );                        BetaHv(i) =1/( exp( (30-V(i))/10 )+1 );

V(i+1)  =V(i) + d* ( Iext - G_Na*m(i)^a*h(i)^b*(V(i)-Ena) - G_K*n(i)^c*(V(i)-Ek) - G_L.*(V(i)-El) );
m(i+1)  =m(i) + dt*( (1-m(i))*AlphaMv(i) - m(i)*BetaMv(i) );
n(i+1)  =n(i) + dt*( (1-n(i))*AlphaNv(i) - n(i)*BetaNv(i) );
h(i+1)  =h(i) + dt*( (1-h(i))*AlphaHv(i) - h(i)*BetaHv(i) );
end
AlphaMv(nt)=( (25-V(nt))/10 )/( exp( (25-V(nt))/10 )-1 );    BetaMv(nt)=4*exp(-V(nt)/18);
AlphaNv(nt)=(10-V(nt))/( 100*( exp( (10-V(nt))/10 )-1) );    BetaNv(nt)=0.125*exp(-V(nt)/80);
AlphaHv(nt)=0.07*exp( -V(nt)/20 );                           BetaHv(nt)=1/( exp( (30-V(nt))/10 )+1 );
