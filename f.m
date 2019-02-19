function fI=f(h)
%Sam Furness
%2/18/19
%Takes concentration values for mRNA or protein and 
    %calculates hill function
n1=2;
n2=2;
n3=2;
K1=0.3*(1000)*10/(3*1.7);%converted to umol/gDW
K2=0.3*(1000)*10/(3*1.7);
K3=0.3*(1000)*10/(3*1.7);
fI=(h^n1)/(K1^n1+h^n1);
end