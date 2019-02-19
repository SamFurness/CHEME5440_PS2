%CHEME5440HW2 Part C
%Sam Furness
%2/17/19

%Define Parameters: Note that rates are converted to 
    %per hour so that the discretization works with 
    %numbers less close to zero
%TRANSCRIPTION
LX1=1200;
LX2=2400;
LX3=600;
LX=1000;
eX=42;
kEX1=3600*eX/LX1;%hr^-1
kEX2=3600*eX/LX2;
kEX3=3600*eX/LX3;
GDW=(1/(6.022*10^23))*(1/1.7)*(10/3)*(1/(6.7*10^-16))*10^6;
%converts molec/cell to umol/gDW
RXT=8000*GDW;%umol/gDW
Gj=200*GDW;%umol/gDW
kmin1=3600*0.1;%hr^-1
kmin2=3600*0.1;
kmin3=3600*0.1;
kplus1=3600*5;%uM^-1hr^-1
kplus2=3600*5;
kplus3=3600*5;
kI1=3600*1/42;%hr^-1
kI2=3600*1/42;
kI3=3600*1/42;
kA1=0;%hr^-1
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
kdX1=3600*log(2)/120;%hr^-1
kdX2=3600*log(2)/120;
kdX3=3600*log(2)/120;
mu=3600*log(2)/(30*60);%hr^-1

%TRANSLATION
LL1=LX1/3;
LL2=LX2/3;
LL3=LX3/3;
LL=333;
eL=14.5;
kEL1=3600*eL/LL1;
kEL2=3600*eL/LL2;
kEL3=3600*eL/LL3;%hr^-1
GDW=(1/(6.022*10^23))*(1/1.7)*(10/3)*(1/(6.7*10^-16))*10^6;
%converts molec/cell to umol/gDW
RLT=20100*GDW;%umol/gDW
kminL1=3600*0.1;%s^-1
kminL2=3600*0.1;
kminL3=3600*0.1;
kplusL1=3600*5;%uM^-1hr^-1
kplusL2=3600*5;
kplusL3=3600*5;
kIL1=3600*1/15;%hr^-1
kIL2=3600*1/15;
kIL3=3600*1/15;
kAL1=0;%hr^-1
kAL2=0;
kAL3=0;
KL1= ((kminL1+kIL1)/kplusL1)*(1000)*10/(3*1.7);
KL2= ((kminL2+kIL2)/kplusL2)*(1000)*10/(3*1.7);
KL3= ((kminL3+kIL3)/kplusL3)*(1000)*10/(3*1.7);
tauL1= (kAL1+kEL1)/kIL1;
tauL2= (kAL2+kEL2)/kIL2;
tauL3= (kAL3+kEL3)/kIL3;
kdL1=3600*log(2)/(20*3600);%hr^-1
kdL2=3600*log(2)/(20*3600);
kdL3=3600*log(2)/(20*3600);

%Define Reaction Rates
%TRANSCRIPTION:
rX1=kEX1*RXT*(Gj/(KX1*tau1+Gj*(tau1+1)));
rX2=kEX2*RXT*(Gj/(KX2*tau2+Gj*(tau2+1)));
rX3=kEX3*RXT*(Gj/(KX3*tau3+Gj*(tau3+1)));

%Setup Differential Equations Matrices
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

%Start by defining a time vector
T=0.01/60;%Time step in minutes
tsim=linspace(0,360,((6/T)+1));
x=zeros(6,length(tsim));%preallocate x vector

%Set initial conditions
x(:,1)=[0;0;0;0;0;0];%Initially no mRNA or protein
m1_0=x(1,1);
m2_0=x(2,1);
m3_0=x(3,1);
p1_0=x(4,1);
p2_0=x(5,1);
p3_0=x(6,1);
I_0=10*(1000)*10/(3*1.7);
%Initial values for r. Uses function f for fI values
u1_0=(WRT1+(WI*f(I_0)))/(1+WRT1+(WI*f(I_0)));
u2_0=(WRT2+(W12*f(p1_0))+(W32*f(p3_0)))/(1+WRT2+(W12*f(p1_0))+(W32*f(p3_0)));
u3_0=(WRT3+(W13*f(p1_0))+(W23*f(p2_0)))/(1+WRT3+(W13*f(p1_0))+(W23*f(p2_0)));
rL1_0=kEL1*RLT*(m1_0/(KL1*tauL1+m1_0*(tauL1+1)));
rL2_0=kEL2*RLT*(m2_0/(KL2*tauL2+m2_0*(tauL2+1)));
rL3_0=kEL3*RLT*(m3_0/(KL3*tauL3+m3_0*(tauL3+1)));
TX1_0=rX1*u1_0;
TX2_0=rX2*u2_0;
TX3_0=rX3*u3_0;
TL4_0=rL1_0;
TL5_0=rL2_0;
TL6_0=rL3_0;

%Set up r vector initial conditions
r0=[TX1_0;TX2_0;TX3_0;TL4_0;TL5_0;TL5_0];
r=zeros(6,length(tsim));
r(1,1)=r0(1);
r(2,1)=r0(2);
r(3,1)=r0(3);
r(4,1)=r0(4);
r(5,1)=r0(5);
r(6,1)=r0(6);

%Performing the Discretization
A_hat=eye(6);
A_hat(1,1)=exp(A(1,1)*T);
A_hat(2,2)=exp(A(2,2)*T);
A_hat(3,3)=exp(A(3,3)*T);
A_hat(4,4)=exp(A(4,4)*T);
A_hat(5,5)=exp(A(5,5)*T);
A_hat(6,6)=exp(A(6,6)*T);
A_inverse=inv(A);
A_sub=A_hat-eye(6);
S_sub=inv(A)*A_sub;
S_hat=S_sub*S;

%Solve Discretization
for i=1:(length(tsim)-1)
   h=tsim(i);
   %Inducer profiles
   if h<60
       Ix=10*(1000)*10/(3*1.7);
   else
       Ix=0*(1000)*10/(3*1.7);
   end
   
   %Get the current r and x vectors
   xk=zeros(6,1);
   rk=zeros(6,1);
   xk(1)=x(1,i);
   xk(2)=x(2,i);
   xk(3)=x(3,i);
   xk(4)=x(4,i);
   xk(5)=x(5,i);
   xk(6)=x(6,i);
  
   rk(1)=r(1,i);
   rk(2)=r(2,i);
   rk(3)=r(3,i);
   rk(4)=r(4,i);
   rk(5)=r(5,i);
   rk(6)=r(6,i);
   
   %Calculate the new entry:
   xkplus1=A_hat*xk+S_hat*rk;
   
   %Update x with new values
   x(1,i+1)=xkplus1(1);
   x(2,i+1)=xkplus1(2);
   x(3,i+1)=xkplus1(3);
   x(4,i+1)=xkplus1(4);
   x(5,i+1)=xkplus1(5);
   x(6,i+1)=xkplus1(6);
   
   %Update r 
   m1=xkplus1(1);
   m2=xkplus1(2);
   m3=xkplus1(3);
   p1=xkplus1(4);
   p2=xkplus1(5);
   p3=xkplus1(6);
   u1=(WRT1+(WI*f(Ix)))/(1+WRT1+(WI*f(Ix)));
   u2=(WRT2+(W12*f(p1))+(W32*f(p3)))/(1+WRT2+(W12*f(p1))+(W32*f(p3)));
   u3=(WRT3+(W13*f(p1))+(W23*f(p2)))/(1+WRT3+(W13*f(p1))+(W23*f(p2)));
   rL1=kEL1*RLT*(m1/(KL1*tauL1+m1_0*(tauL1+1)));
   rL2=kEL2*RLT*(m2/(KL2*tauL2+m2_0*(tauL2+1)));
   rL3=kEL3*RLT*(m3/(KL3*tauL3+m3_0*(tauL3+1)));
   TX1=rX1*u1;
   TX2=rX2*u2;
   TX3=rX3*u3;
   TL4=rL1;
   TL5=rL2;
   TL6=rL3;
%now r
   r(1,i+1)=TX1;
   r(2,i+1)=TX2;
   r(3,i+1)=TX3;
   r(4,i+1)=TL4;
   r(5,i+1)=TL5;
   r(6,i+1)=TL6;
end
      
figure(1)
I_plot=zeros(1,length(tsim));
I_plot(1:6000)=10*50;%scale appropriately
hold on
plot(tsim,x(4,:));
plot(tsim,x(5,:));
plot(tsim,x(6,:));
plot(tsim,I_plot); % comment this out if you don't want inducer
title('Memory Circuit Proteins Over Time-Discretization Method')
xlabel('Time (min)')
ylabel('Concentration (umol/gDW)')
legend('P1','P2','P3','Inducer (scaled)','P1-Discretized','P2-Discretized','P3-Discretized')

hold off
    
    
    
    
    
    