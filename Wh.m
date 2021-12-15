function y = Wh(bf,bm,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im)
%Wp calculates the fitness of the host at equilibrium Sf, Sm, If, Im
%   gammas and bf and bm (formerly just b) are _mutant_ values 
F=zeros(4);
Vinv=zeros(4);
N=Sf+Sm+If+Im;

% populate F column by column
F(1,1)= bf * (Sm + Im) / N;
F(2,1)= bf * (Sm + Im) / N;

F(1,2)= bm * (Sf + (1-v)*If) / N;
F(2,2)= bm * (Sf + (1-v)*If) / N;
F(3,2)= bm * v * If / N;
F(4,2)= bm * v * If / N;

F(1,3)= (1-v) * bf * (Sm + Im) / N;
F(2,3)= (1-v) * bf * (Sm + Im) / N;
F(3,3)= v * bf * (Sm + Im) / N;
F(4,3)= v * bf * (Sm + Im) / N;

F(1,4)= bm * (Sf + (1-v)*If) / N;
F(2,4)= bm * (Sf + (1-v)*If) / N;
F(3,4)= bm * v * If / N;
F(4,4)= bm * v * If / N;

F = 0.5 * F;

% populate Vinv column by column
tauf=( mu*N*gammaf + (mu*N+alphaf) * (mu*N + betaff*If + betafm*Im) )^(-1);
taum=( mu*N*gammam + (mu*N+alpham) * (mu*N + betamf*If + betamm*Im) )^(-1);

Vinv(1,1)= ( mu * N + alphaf + gammaf ) * tauf;
Vinv(3,1)= ( betaff * If + betafm * Im ) * tauf;

Vinv(2,2)= ( mu * N + alpham + gammam ) * taum;
Vinv(4,2)= ( betamf * If + betamm * Im ) * taum;

Vinv(1,3)= gammaf * tauf;
Vinv(3,3)= ( mu * N + betaff * If + betafm * Im ) * tauf;

Vinv(2,4)= gammam * taum;
Vinv(4,4)= ( mu * N + betamf * If + betamm * Im ) * taum;

y = max(abs(eig(F*Vinv)));
end

