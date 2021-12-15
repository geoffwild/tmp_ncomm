% demo that finds converge-stable and ESS values of gammaf and gammaf
% ResEquil.m and Wh.m should be in the same directory
clear all

% set parameter values 
mu=0.1;
v=0.2;
alphaf = 0.95; 
alpham = 1.05; 

% establish tradeoff faced by pathogen
betaXf =@(a) 1.1 * a./(a + 0.5);
betaXm =@(a) 1.1 * a./(a + 0.5);

betaff=betaXf(alphaf);
betamf=betaXf(alphaf);
    
betafm=betaXm(alpham);
betamm=betaXm(alpham);

% establish tradeoff faced by host
pf =@(gamma) exp( -0.1 * gamma.^2 );
pm =@(gamma) exp( -0.1 * gamma.^2 );


% begin with a poor guess at ESS value
gammaf = rand();
gammam = rand();

% update guesses iteratively to find convergence-stable pair of gammas
tol = 1e-06;
maxiter = 1e05;
iter=0;

grad=ones(2,1);
while norm(grad) > tol
    iter = iter + 1;
    if iter > maxiter
        break
    end
    
    % establish fecundity b
    bmax = 2;
    b = bmax * pf(gammaf) * pm(gammam);
    
    % determine resident equilibrium values
    [y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
    Sf = y(1);
    Sm = y(2);
    If = y(3);
    Im = y(4);
    
    % determine mutant traits
    delta = 1e-01;
    gfup = gammaf + 0.5 * delta;
    gfdn = gammaf - 0.5 * delta;
    gmup = gammam + 0.5 * delta;
    gmdn = gammam - 0.5 * delta;

    while or(gfdn<0, gmdn<0) 
        delta = 0.1 * delta;
        gfup = gammaf + 0.5 * delta;
        gfdn = gammaf - 0.5 * delta;
        gmup = gammam + 0.5 * delta;
        gmdn = gammam - 0.5 * delta;
    end

    bfup = bmax * pf(gfup) * pm(gammam);
    bfdn = bmax * pf(gfdn) * pm(gammam);
 
    bmup = bmax * pf(gammaf) * pm(gmup);
    bmdn = bmax * pf(gammaf) * pm(gmdn);
    
    % determine selection gradient acting on gammaf
    Whfup = Wh(bfup,b,mu,v,gfup,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
    Whfdn = Wh(bfdn,b,mu,v,gfdn,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
    grad(1) = Whfup - Whfdn;
    
    % determine selection gradient acting on gammam
    Whmup = Wh(b,bmup,mu,v,gammaf,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
    Whmdn = Wh(b,bmdn,mu,v,gammaf,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
    grad(2) = Whmup - Whmdn;
    
    step = 1;
    gammaf = gammaf + step * grad(1);
    gammam = gammam + step * grad(2);
end

display(gammaf);
display(gammam);
display(iter);

% check ESS

% establish fecundity(redundant here but can cut paste)
b = bmax * pf(gammaf) * pm(gammam);

% determine resident equilibrium values (redundant here but can cut paste)
[y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
Sf = y(1);
Sm = y(2);
If = y(3);
Im = y(4);
    
% determine mutant traits
delta = 1e-01;
gfup = gammaf + 0.5 * delta;
gfdn = gammaf - 0.5 * delta;
gmup = gammam + 0.5 * delta;
gmdn = gammam - 0.5 * delta;

while or(gfdn<0, gmdn<0) 
    delta = 0.1 * delta;
    gfup = gammaf + 0.5 * delta;
    gfdn = gammaf - 0.5 * delta;
    gmup = gammam + 0.5 * delta;
    gmdn = gammam - 0.5 * delta;
end

bfup = bmax * pf(gfup) * pm(gammam);
bfdn = bmax * pf(gfdn) * pm(gammam);

bmup = bmax * pf(gammaf) * pm(gmup);
bmdn = bmax * pf(gammaf) * pm(gmdn);

Wh0 = Wh(b,b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
% Wp0 should be 1, which can be verified by un-commenting next line
display(Wh0);

Whfup = Wh(bfup,b,mu,v,gfup,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whfdn = Wh(bfdn,b,mu,v,gfdn,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whmup = Wh(b,bmup,mu,v,gammaf,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whmdn = Wh(b,bmdn,mu,v,gammaf,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);

Whupup= Wh(bfup,bmup,mu,v,gfup,gmup,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
Whdndn= Wh(bfdn,bmdn,mu,v,gfdn,gmdn,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);

H=zeros(2);
H(1,1) = (Whfup - 2 * Wh0 + Whfdn) / ((0.5 * delta)^2);
H(2,2) = (Whmup - 2 * Wh0 + Whmdn) / ((0.5 * delta)^2);
H(1,2) = (Whupup - Whfup - Whmup + 2 * Wh0 - Whfdn - Whmdn + Whdndn) / ( 2 * 0.5 * delta * 0.5 *delta );
H(2,1) = (Whupup - Whfup - Whmup + 2 * Wh0 - Whfdn - Whmdn + Whdndn) / ( 2 * 0.5 * delta * 0.5 *delta ); 
evals = eig(H);
if and(evals(1)<0, evals(2)<0)
    display('ESS');
else
    display('not ESS');
end
