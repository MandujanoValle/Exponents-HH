%%%%%%%%             To understand the code read article        %%%%%%%%%%%%
%%  PARAMETER IDENTIFICATION PROBLEM IN THE HODGKIN AND HUXLEY MODEL      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    ESTIMATE EXPONENTS a, b, c                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Model of Hodgkin and Huxley Temporal                    %
%                Euler Explicit:  To estimate  a, b, c                     %
%                                                                          %
%|--------------------------------                                         %
%| CV_t=Iext-G_Na.m^a.h^b.(V-E_Na)-G_K.n^c.(V-E_K)-G_L.(V-E_L)             %
%| m_t=(1-m).AlphaNm(V)-m.BetaNm(V)                                        %
%| h_t=(1-h).AlphaNh(V)-h.BetaNh(V)                                        %
%| n_t=(1-n).AlphaNn(V)-n.BetaNn(V)                                        %
%| V(0)=V_0;   m(0)=m_0    h(0)=h_0   n(0)=n_0                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clc;close all;clear all 
global G_Na G_K G_L Ena Ek El nt dt  m n h V Vp  C_M Iext                    %|
%%-- Start: defining Ordinari Diferential Equation (ODE (7) of paper)   --%%
                                                                           %|
%----    Set the time interval  ([0 T] in the paper) ti=0 and tf=T      ---%|
ti=0;     tf=5;   nt=500;                                                %|
t=linspace(ti,tf,nt) ;  dt = t(2)-t(1) ;                                   %|
                                                                           %|
%----                  known parameters                                 ---%|
C_M=1; Iext=0; Ena=115; Ek=-12; El=10.598; G_Na=120; G_K=36; G_L=0.3;      %|
                                                                           %|
%----                  initial conditions                               ---%|
V=zeros(1,nt); m=zeros(1,nt); n=zeros(1,nt); h=zeros(1,nt);                %|
V(1)=-25;        m(1)=0.5;      n(1)=0.4;      h(1)=0.4;                   %|
                                                                           %|
%---                unknown parameters ( Exponents )                       %|
a =3;          b =1;            c =4;                                      %|
                                                                           %|
%---                guess initial ( Exponents )                            %|
ak=0;          bk=0;            ck=0;                                      %|
                                                                           %|
%---   perturbation of the electrical potential (in percentage )        ---%|
MaxErro=1/100;                                                             %|
                                                                           %|
%---             for the stop criterion                                 ---%|
tau=1.01;                                                                  %|
                                                                           %|
%%---                              End                                ---%%%|

%-----------            Calculating the exact Vexa              ------------%
Vexa=Vexata(a,b,c);

%----------        Making the pertubation of Vexa in Vp         ------------%
Vp=Vexa + (-MaxErro+(2*MaxErro).*rand(1,nt)).*Vexa;

%---------  Calculing delta for the equation (9) of the paper   ------------%
delta=MaxErro*sqrt( dt^1*sum( (Vexa).^2 ) );


%%--------------------       k=========1       ------------%%
k=0;

while(0==k || tau*delta<=ResiduoV(k)) 
%for  i=1:10000
k=k+1;  
%-------------           Calculing  Vk, Uk, mk, nk, hk               ------------%
  [Vk,Uk,mk,nk,hk]=Iaprox(ak,bk,ck);

%-------------          Calculing of the residue: ||Vp-Vk||         -------------%
  ResiduoV(k)=sqrt( dt*sum( (Vp-Vk).^2 ) );

%-------------               Calculing the Error                  ---------------%
  Error(k)=sqrt( [(ak-a)^2+(bk-b)^2+(ck-c)^2 ]/[(a)^2+(b)^2+(c)^2] )*100;
  aa(k)=ak;    bb(k)=bk;    cc(k)=ck;  

%------                  To the minimal error method           --------------%
  Wk=sum( dt*(Vp-Vk).^2 )/...
        ( ( sum( G_Na.*(Vk-Ena).*mk.^ak.*hk.^bk.*Uk.*log(mk) ) )^2+...
          ( sum( G_Na.*(Vk-Ena).*mk.^ak.*hk.^bk.*Uk.*log(hk) ) )^2+...
          ( sum( G_K .*(Vk-Ek ).*nk.^ck        .*Uk.*log(nk) ) )^2  );
%  Wk=1;
%-------------                      Print                         --------------%
fprintf('%10.5f\t',k,ak,bk,ck,Error(k),ResiduoV(k),Wk);fprintf('\n\n\n');

%------------         Calculing the iteration k+1             --------------%
ak=ak+Wk*dt*sum( G_Na.*(Vk-Ena).*mk.^ak.*hk.^bk.*Uk.*log(mk) ); %For a
bk=bk+Wk*dt*sum( G_Na.*(Vk-Ena).*mk.^ak.*hk.^bk.*Uk.*log(hk) ); %For b
ck=ck+Wk*dt*sum( G_K .*(Vk-Ek ).*nk.^ck        .*Uk.*log(nk) ); %For c
end

k=linspace(1,k,k);  ks=k(1:100:max(k)); if (max(k)~= max(ks) ) k= [ks max(k) ]; end;
ak=aa(k); bk=bb(k);  ck=cc(k); ResiduoV= ResiduoV(k); Error=Error(k);
save('Example-11.txt','Vexa','Vp','t','-ascii');
save('Example-12.txt','ak','bk','ck','ResiduoV','Error','k','-ascii');
