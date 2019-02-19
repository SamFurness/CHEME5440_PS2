%Problem 2 on PS2 of CHEME 5440
%Sam Furness
%2/12/19
function dxdt = CHEME5440HW2B(t,x,Ix)
%Define Parameters:
m1=x(1);
m2=x(2);
m3=x(3);
p1=x(4);
p2=x(5);
p3=x(6);
%TRANSCRIPTION
LX1=1200;
LX2=2400;
LX3=600;
LX=1000;
eX=42;
kEX1=eX/LX1;
kEX2=eX/LX2;
kEX3=eX/LX3;%s^-1
GDW=(1/(6.022*10^23))*(1/1.7)*(10/3)*(1/(6.7*10^-16))*10^6;
%converts molec/cell to umol/gDW
RXT=8000*GDW;%umol/gDW
Gj=200*GDW;%umol/gDW
kmin1=0.1;%s^-1
kmin2=0.1;
kmin3=0.1;
kplus1=5;%uM^-1s^-1
kplus2=5;
kplus3=5;
kI1=1/42;%s^-1
kI2=1/42;
kI3=1/42;
kA1=0;%s^-1
kA2=0;
kA3=0;
KX1= ((kmin1+kI1)/kplus1)*(1000)*10/(3*1.7);
KX2= ((kmin2+kI2)/kplus2)*(1000)*10/(3*1.7);
KX3= ((kmin3+kI3)/kplus3)*(1000)*10/(3*1.7);
tau1= (kA1+kEX1)/kI1;
tau2= (kA2+kEX2)/kI2;
tau3= (kA3+kEX3)/kI3;
WRT1=0.00001;
WI=100;
WRT2=0.00001;
W12=100;
W32=10;
WRT3=0.00001;
W13=100;
W23=10;
n1=2;
n2=2;
n3=2;
K1=0.3*(1000)*10/(3*1.7);%converted to umol/gDW
K2=0.3*(1000)*10/(3*1.7);
K3=0.3*(1000)*10/(3*1.7);
kdX1=log(2)/120;%s^-1
kdX2=log(2)/120;
kdX3=log(2)/120;
mu=log(2)/(30*60);%s^-1

%TRANSLATION
LL1=LX1/3;
LL2=LX2/3;
LL3=LX3/3;
LL=333;
eL=14.5;
kEL1=eL/LL1;
kEL2=eL/LL2;
kEL3=eL/LL3;%s^-1
GDW=(1/(6.022*10^23))*(1/1.7)*(10/3)*(1/(6.7*10^-16))*10^6;
%converts molec/cell to umol/gDW
RLT=20100*GDW;%umol/gDW
kminL1=0.1;%s^-1
kminL2=0.1;
kminL3=0.1;
kplusL1=5;%uM^-1s^-1
kplusL2=5;
kplusL3=5;
kIL1=1/15;%s^-1
kIL2=1/15;
kIL3=1/15;
kAL1=0;%s^-1
kAL2=0;
kAL3=0;
KL1= ((kminL1+kIL1)/kplusL1)*(1000)*10/(3*1.7);
KL2= ((kminL2+kIL2)/kplusL2)*(1000)*10/(3*1.7);
KL3= ((kminL3+kIL3)/kplusL3)*(1000)*10/(3*1.7);
tauL1= (kAL1+kEL1)/kIL1;
tauL2= (kAL2+kEL2)/kIL2;
tauL3= (kAL3+kEL3)/kIL3;
kdL1=log(2)/(20*3600);%s^-1
kdL2=log(2)/(20*3600);
kdL3=log(2)/(20*3600);

%Define Control Functions and Inducer Vector
I=Ix*(1000)*10/(3*1.7);%convert to umol/gDW
fI=(I^n1)/(K1^n1+I^n1);
fI12=(p1^n2)/(K2^n2+p1^n2);
fI32=(p3^n2)/(K2^n2+p3^n2);
fI13=(p1^n3)/(K3^n3+p1^n3);
fI23=(p2^n3)/(K3^n3+p2^n3);
u1=(WRT1+(WI*fI))/(1+WRT1+(WI*fI));
u2=(WRT2+(W12*fI12)+(W32*fI32))/(1+WRT2+(W12*fI12)+(W32*fI32));
u3=(WRT3+(W13*fI13)+(W23*fI23))/(1+WRT3+(W13*fI13)+(W23*fI23));

%Define Reaction Rates
%TRANSCRIPTION:
rX1=kEX1*RXT*(Gj/(KX1*tau1+Gj*(tau1+1)));
rX2=kEX2*RXT*(Gj/(KX2*tau2+Gj*(tau2+1)));
rX3=kEX3*RXT*(Gj/(KX3*tau3+Gj*(tau3+1)));
%TRANSLATION:
rL1=kEL1*RLT*(m1/(KL1*tauL1+m1*(tauL1+1)));
rL2=kEL2*RLT*(m2/(KL2*tauL2+m2*(tauL2+1)));
rL3=kEL3*RLT*(m3/(KL3*tauL3+m3*(tauL3+1)));

%Setup Differential Equations Matrices
vec=[m1;m2;m3;p1;p2;p3];
r=[rX1*u1;rX2*u2;rX3*u3;rL1;rL2;rL3];
A=[-(mu+kdX1),0,0,0,0,0;
   0,-(mu+kdX2),0,0,0,0;
   0,0,-(mu+kdX3),0,0,0;
   0,0,0,-(mu+kdL1),0,0;
   0,0,0,0,-(mu+kdL2),0;
   0,0,0,0,0,-(mu+kdL3)];
S=[1,0,0,0,0,0;
   0,1,0,0,0,0;
   0,0,1,0,0,0;
   0,0,0,1,0,0;
   0,0,0,0,1,0;
   0,0,0,0,0,1];
%Create differential equation:
dxdt = (A*vec)+(S*r);

end










