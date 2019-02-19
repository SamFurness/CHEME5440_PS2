%Running file for CHEME5440 HW2. ODE

%Set up time vector
tStart = 0.0;
tStop = 3600.0;
tSim = [tStart,tStop];
%Set up initial conditions
x0 = [0.0; %m1
      0.0; %m2
      0.0; %m3
      0.0; %p1
      0.0; %p2
      0.0];%p3
I1=10;
%Run ode for while inducer present:
[t,X]=ode45(@(t,x) CHEME5440HW2B(t,x,I1),tSim,x0);

%Set up new time vector
tStart2 = 0;
tStop2 = 300*60;
tSim2 = [tStart2,tStop2];
%Set up new initial conditions
x02=[X(end,1);
     X(end,2);
     X(end,3);
     X(end,4);
     X(end,5);
     X(end,6)];
I2=0;
[t2,Y]=ode45(@(t2,x2) CHEME5440HW2B(t2,x2,I2),tSim2,x02);

%Concatenate two vectors:
t3=t(end)+t2;
t_tot=cat(1,t,t3);
X_tot=cat(1,X,Y);

figure(1)
I_plot=zeros(1,length(t_tot));
I_plot(1:length(t))=I1*50;
figure(1)
hold on
t_plot=(1/60)*t_tot;
plot(t_plot,X_tot(:,4))
plot(t_plot,X_tot(:,5))
plot(t_plot,X_tot(:,6))
plot(t_plot,I_plot)
title('Memory Circuit Proteins Over Time')
xlabel('Time (min)')
ylabel('Concentration (umol/gDW)')
legend('P1','P2','P3','Inducer (scaled)')
hold off



