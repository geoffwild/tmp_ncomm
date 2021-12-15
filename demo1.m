% population dynamics demo
% ResRzero.m and ResEquil.m should be in the same directory
clear all

% set parameter values 
% (sex specific parameters have been equated)
b=2;
mu=0.1;
v=0.2;
gammaf=1.5;
gammam=1.5;
alphaf=1;
alpham=1;
betaff=0.9;
betafm=0.9;
betamf=0.9;
betamm=0.9;

% calculate the equilibrium state and the vector of derivatives in order
% to verify that derivatives are within a tolerance of zero
[y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
display(y)
display(dydt)

% calculate R0 and verify that it is greater than 1 when infections 
% occur at non-zero densities at equilibrium
R0=ResRzero(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
display(R0)

% create a plot to illustrate the transcritical bifurcation
% at R0 = 1
nobs = 1000;
beta = linspace(0,1,nobs);
hData = nan(nobs,1);
vSfData = nan(nobs,1);
vSmData = nan(nobs,1);
vIfData = nan(nobs,1);
vImData = nan(nobs,1);
for i=1:nobs
    betaff=beta(i);
    betafm=beta(i);
    betamf=beta(i);
    betamm=beta(i);
    [y,dydt]=ResEquil(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
    vSfData(i) = y(1); 
    vSmData(i) = y(2);
    vIfData(i) = y(3);
    vImData(i) = y(4);
    hData(i)=ResRzero(b,mu,v,gammaf,gammam,alphaf,alpham,betaff,betafm,betamf,betamm);
end

figure 
hold on
plot(hData,vSfData, '-m', 'LineWidth', 2)
plot(hData,vSmData, '--b', 'LineWidth', 1)
plot(hData,vIfData, '-c', 'LineWidth', 2)
plot(hData,vImData, '--r', 'LineWidth', 1)
xlabel('R_0')
ylabel('equilibrium numbers of individuals')
legend('S_f', 'S_m', 'I_f', 'I_m')