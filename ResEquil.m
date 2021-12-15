function [y,dydt] = ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm)
%ResEquil finds the equilibrium of the resident system of ODEs iteratively
%   Detailed explanation goes here

dydt = ones(4,1);
y= [0.25*b/mu; 0.25*b/mu; 0.1*rand(); 0.1*rand()];
step=1e-01;
tol=1e-06;

while norm(dydt) > tol
    Sf = y(1);
    Sm = y(2);
    If = y(3);
    Im = y(4);
    N = Sf+Sm+If+Im;

    dydt(1) = b*(Sf + (1-v)*If)*(Sm+Im)/N + gammaf*If - mu*N*Sf - Sf*betaff*If - Sf*betafm*Im; 
    dydt(2) = b*(Sf + (1-v)*If)*(Sm+Im)/N + gammam*Im - mu*N*Sm - Sm*betamf*If - Sm*betamm*Im;
    dydt(3) = b*v*If*(Sm+Im)/N - gammaf*If - mu*N*If - alphaf*If + Sf*betaff*If + Sf*betafm*Im; 
    dydt(4) = b*v*If*(Sm+Im)/N - gammam*Im - mu*N*Im - alpham*Im + Sm*betamf*If + Sm*betamm*Im; 
    
    y=y+step*dydt;
end
end

