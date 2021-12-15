function y = Wp(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im)
%Wp calculates the fitness of the parasite at equilibrium Sf, Sm, If, Im
%   alphas and betas are _mutant_ values 
K = zeros(2);
N=Sf+Sm+If+Im;
K(1,1) = (b*v*(Sm+Im)/N + Sf*betaff) / (N*mu + gammaf + alphaf);
K(2,1) = (b*v*(Sm+Im)/N + Sm*betamf) / (N*mu + gammaf  +alphaf);
K(1,2) = Sf*betafm / (N*mu + gammam + alpham);
K(2,2) = Sm*betamm / (N*mu + gammam + alpham);
y = max(abs(eig(K)));
end

