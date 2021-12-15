% demo that finds converge-stable and ESS values of alphaf and alphaf
% ResEquil.m and Wp.m should be in the same directory
clear all

% set parameter values 
% (sex specific parameters have been equated)
b=2;
mu=0.1;
v=0.2;
gammaf=1.5;
gammam=1.5;

% establish tradeoff faced by pathogen
% (still the same for sexes)
betaXf =@(a) 1.1 * a./(a + 0.5);
betaXm =@(a) 1.1 * a./(a + 0.5);

% begin with a poor guess at ESS value
alphaf = rand();
alpham = rand();

% update guesses iteratively to find convergence-stable pair of alphas
tol = 1e-06;
maxiter = 1e05;
iter=0;

grad=ones(2,1);
while norm(grad) > tol
    iter = iter + 1;
    if iter > maxiter
        break
    end
    
    % establish transmissibilities
    betaff=betaXf(alphaf);
    betamf=betaXf(alphaf);
    
    betafm=betaXm(alpham);
    betamm=betaXm(alpham);
    
    % determine resident equilibrium values
    [y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
    Sf = y(1);
    Sm = y(2);
    If = y(3);
    Im = y(4);
    
    % determine mutant traits
    delta = 1e-01;
    afup = alphaf + 0.5 * delta;
    afdn = alphaf - 0.5 * delta;
    amup = alpham + 0.5 * delta;
    amdn = alpham - 0.5 * delta;

    while or(afdn<0, amdn<0) 
        delta = 0.1 * delta;
        afup = alphaf + 0.5 * delta;
        afdn = alphaf - 0.5 * delta;
        amup = alpham + 0.5 * delta;
        amdn = alpham - 0.5 * delta;
    end

    betaffup=betaXf(afup);
    betaffdn=betaXf(afdn); 

    betamfup=betaXf(afup);
    betamfdn=betaXf(afdn);

    betafmup=betaXm(amup);
    betafmdn=betaXm(amdn);

    betammup=betaXm(amup);
    betammdn=betaXm(amdn);
    
    % determine selection gradient acting on alphaf
    Wpfup = Wp(b,mu,v,gammaf,gammam,afup,alpham,betaffup,betafm,betamfup,betamm,Sf,Sm,If,Im);
    Wpfdn = Wp(b,mu,v,gammaf,gammam,afdn,alpham,betaffdn,betafm,betamfdn,betamm,Sf,Sm,If,Im);
    grad(1) = Wpfup - Wpfdn;
    
    % determine selection gradient acting on alpham
    Wpmup = Wp(b,mu,v,gammaf,gammam,alphaf,amup,betaff,betafmup,betamf,betammup,Sf,Sm,If,Im);
    Wpmdn = Wp(b,mu,v,gammaf,gammam,alphaf,amdn,betaff,betafmdn,betamf,betammdn,Sf,Sm,If,Im);
    grad(2) = Wpmup - Wpmdn;
    
    step = 1;
    alphaf = alphaf + step * grad(1);
    alpham = alpham + step * grad(2);
end

% as expected, alphaf tends to a lower value than alpham
display(alphaf);
display(alpham);
display(iter);

% check ESS

% establish transmissibilities (redundant here but can cut paste)
betaff=betaXf(alphaf);
betamf=betaXf(alphaf);

betafm=betaXm(alpham);
betamm=betaXm(alpham);

% determine resident equilibrium values (redundant here but can cut paste)
[y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
Sf = y(1);
Sm = y(2);
If = y(3);
Im = y(4);
    
% determine mutant traits
delta = 1e-01;
afup = alphaf + 0.5 * delta;
afdn = alphaf - 0.5 * delta;
amup = alpham + 0.5 * delta;
amdn = alpham - 0.5 * delta;

while or(afdn<0, amdn<0) 
    delta = 0.1 * delta;
    afup = alphaf + 0.5 * delta;
    afdn = alphaf - 0.5 * delta;
    amup = alpham + 0.5 * delta;
    amdn = alpham - 0.5 * delta;
end

betaffup=betaXf(afup);
betaffdn=betaXf(afdn); 

betamfup=betaXf(afup);
betamfdn=betaXf(afdn);

betafmup=betaXm(amup);
betafmdn=betaXm(amdn);

betammup=betaXm(amup);
betammdn=betaXm(amdn);

Wp0 =  Wp(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm,Sf,Sm,If,Im);
% Wp0 should be 1, which can be verified by un-commenting next line
% display(Wp0);

Wpfup = Wp(b,mu,v,gammaf,gammam,afup,alpham,betaffup,betafm,betamfup,betamm,Sf,Sm,If,Im);
Wpfdn = Wp(b,mu,v,gammaf,gammam,afdn,alpham,betaffdn,betafm,betamfdn,betamm,Sf,Sm,If,Im);
Wpmup = Wp(b,mu,v,gammaf,gammam,alphaf,amup,betaff,betafmup,betamf,betammup,Sf,Sm,If,Im);
Wpmdn = Wp(b,mu,v,gammaf,gammam,alphaf,amdn,betaff,betafmdn,betamf,betammdn,Sf,Sm,If,Im);

Wpupup = Wp(b,mu,v,gammaf,gammam,afup,amup,betaffup,betafmup,betamfup,betammup,Sf,Sm,If,Im);
Wpdndn = Wp(b,mu,v,gammaf,gammam,afdn,amdn,betaffdn,betafmdn,betamfdn,betammdn,Sf,Sm,If,Im);

H=zeros(2);
H(1,1) = (Wpfup - 2 * Wp0 + Wpfdn) / ((0.5 * delta)^2);
H(2,2) = (Wpmup - 2 * Wp0 + Wpmdn) / ((0.5 * delta)^2);
H(1,2) = (Wpupup - Wpfup - Wpmup + 2 * Wp0 - Wpfdn - Wpmdn + Wpdndn) / ( 2 * 0.5 * delta * 0.5 *delta );
H(2,1) = (Wpupup - Wpfup - Wpmup + 2 * Wp0 - Wpfdn - Wpmdn + Wpdndn) / ( 2 * 0.5 * delta * 0.5 *delta ); 
evals = eig(H);
if and(evals(1)<0, evals(2)<0)
    display('ESS');
else
    display('not ESS');
end