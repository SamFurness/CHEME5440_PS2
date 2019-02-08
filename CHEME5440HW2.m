%Using HW1 to bring up differing plots of mRNA versus I
%Sam Furness
%2/7/19
q=input('What do you want to plot?');
a=0.8;
b=1.2;
kEj=0.0137;%s^-1
kEj_array=[a*kEj,kEj,b*kEj];
RXT=19.8;%uM
RXT_array=[a*RXT,RXT,b*RXT];
Gj=6.196;%uM
Gj_array=[a*Gj,Gj,b*Gj];
kmin=0.1;%s^-1
kmin_array=[a*kmin,kmin,b*kmin];
kplus=5;%uM^-1s^-1
kplus_array=[a*kplus,kplus,kplus*b];
kI=1/42;%s^-1
kI_array=[a*kI,kI,b*kI];
kA=0;%s^-1
kA_array=[a*kA,kA,b*kA];
KXj= (kmin+kI)/kplus;
KXj_array=[a*KXj,KXj,b*KXj];
tau= (kA+kEj)/kI;
tau_array=[a*tau,tau,b*tau];
W1=0.26;
W1_array=[a*W1,W1,b*W1];
W2=300;
W2_array=[a*W2,W2,b*W2];
n=1.5;
n_array=[a*n,n,b*n];
K=0.3;
K_array=[a*K,K,b*K];
kd=log(2)/120;%s^-1
kd_array=[a*kd,kd,b*kd];
mu=1.14/3600;%s^-1
mu_array=[a*mu,mu,b*mu];

I=linspace(0.0001,10,100000);%mM
for i=1:length(I)
    fI(i)=(I(i)^n)/(K^n+I(i)^n);
end
for i=1:length(I)
    for j=1:3
    fI_array(j,i)=(I(i)^n_array(j))/(K_array(j)^n_array(j)+I(i)^n_array(j));
    end
end

for i=1:length(I)
    mj(i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
end

if q=='kEj'
tau= (kA+kEj_array)/kI;
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj_array(j)*RXT*Gj)/(KXj*tau(j)+Gj*(tau(j)+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered kEj')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='RXT';

for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT_array(j)*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer ConcentrationWith Altered RXT')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='Gj_'
    
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj_array(j))/(KXj*tau+Gj_array(j)*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer ConcentrationWith Altered Gj')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')


elseif q=='kmi'
    KXj=(kmin_array+kI)/kplus;
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj(j)*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer ConcentrationWith Altered kmin')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='kpl'
    KXj=(kmin+kI)./kplus_array;
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj(j)*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer ConcentrationWith Altered kplus')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='kI_'
    KXj=(kmin+kI_array)/kplus;
    tau= (kA+kEj)./kI_array;
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj(j)*tau(j)+Gj*(tau(j)+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered kI')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='KXj'
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj_array(j)*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered KXj')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='tau'
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau_array(j)+Gj*(tau_array(j)+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered tau')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='W1_'
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1_array(j)+W2*fI(i))/(1+W1_array(j)+W2*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered W1')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='W2_'
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2_array(j)*fI(i))/(1+W1+W2_array(j)*fI(i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered W2')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='n__'
for i=1:length(I)
    for j=1:3
    fI_array(j,i)=(I(i)^n_array(j))/(K^n_array(j)+I(i)^n_array(j));
    end
end
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI_array(j,i))/(1+W1+W2*fI_array(j,i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered n')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='K__'
for i=1:length(I)
    for j=1:3
    fI_array(j,i)=(I(i)^n)/(K_array(j)+I(i)^n);
    end
end
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI_array(j,i))/(1+W1+W2*fI_array(j,i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered K')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='fI_'
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI_array(j,i))/(1+W1+W2*fI_array(j,i))))/(kd+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered fI')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='kd_'
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd_array(j)+mu);
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered kd')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')

elseif q=='mu_'
for i=1:length(I)
    for j=1:3
    mj_array(j,i)=(((kEj*RXT*Gj)/(KXj*tau+Gj*(tau+1)))*((W1+W2*fI(i))/(1+W1+W2*fI(i))))/(kd+mu_array(j));
    end
end
semilogx(I,mj_array(1,:),I,mj_array(2,:),I,mj_array(3,:))
title('mRNA versus Inducer Concentration With Altered mu')
xlabel('Inducer Concentration mM')
ylabel('mRNA Concentration uM')
legend({'Decrease20%','Original','Increase20%'},'Location','southeast')


end



%logI=log10(I);
%logmj=log10(mj);
%semilogx(I,mj)
%plot(logI,mj);
%plot(I,logmj);
%title('mRNA versus Inducer Concentration')
%xlabel('Inducer Concentration mM')
%ylabel('mRNA Concentration uM')
clear

