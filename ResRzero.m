function R0 = ResRzero(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

K=zeros(2);

S0=0.25*b/mu;

K(1,1) = (b*v/2 + S0*betaff)/(2*S0*mu+alphaf+gammaf);
K(2,1) = (b*v/2 + S0*betamf)/(2*S0*mu+alphaf+gammaf);

K(1,2) = S0*betafm/(2*S0*mu+alpham+gammam);
K(2,2) = S0*betamm/(2*S0*mu+alpham+gammam);

R0 = max(abs(eig(K)));

end

