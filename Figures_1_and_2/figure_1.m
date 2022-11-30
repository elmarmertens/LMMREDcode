%% schwanNK example with inflation and output observed

%% clear workspace
clc
clear
close all

%% load toolboxes
path(pathdef)

addpath ../matlabbox/emtools/
addpath ../matlabbox/emtexbox/
addpath ../matlabbox/emgibbsbox/
addpath ../matlabbox/emeconometrics/
addpath ../matlabbox/emstatespace/

rng(100)
tic
%#ok<*UNRCH>
%#ok<*NOPTS>

%% calibration
% textbook values from Gali:
beta    = .99;
phi     = 1;   % inverse Frisch labor elasticity
sigma   = 1;   % Risk-Aversion / Inverse IES

% Philipps Curve
theta   = .75;    % Calvo Probability of *not* resetting price
kappa   = (1 - theta) * (1 - beta * theta) / theta * (sigma + phi); % lambda  = (1 - theta) * (1 - beta * theta) / theta;

% policy rule:





phi_pi = 2.5;
phi_x  = 0.5;
doTaylorWithNaturalRate = false;  
% choice for doTaylorWithNaturalRate:
% if false: i(t) = phi_pi * pi(t) + phi_x * x(t) 
% if true:  i(t) = rbar(t) + phi_pi * pi(t) + phi_x * x(t) 

% stochastic processes
sig_u  = 0.3;
rho_u  = 0.75;

sig_y  = 0.3;
rho_y  = 0.75;




Ny      = 2;

%% pi-y full info system
% A y(t+1) = B y(t) + D v(t)
% y(t) = [pi, y]

A = [beta 0; 1 sigma];
B = [1, -kappa; phi_pi, phi_x+sigma];
if doTaylorWithNaturalRate
    D = [-1, kappa; 0, -phi_x + sigma * (rho_y - 1)];
else
    D = [-1 kappa; 0 -phi_x];
end

F = [rho_u 0; 0 rho_y];
I = eye(Ny);

Gvec = -(kron(I,B) - kron(F',A)) \ D(:);
G    = reshape(Gvec, Ny, Ny);


checkdiff(A * G * F, B * G + D);

g_piu = G(1,1);
g_piy = G(1,2);
g_yu  = G(2,1);
g_yy  = G(2,2);

%% check coefficients
checkdiff(g_piu, (1 + kappa * g_yu) / (1 - beta * rho_u));
checkdiff(g_yu, - (phi_pi - rho_u) / (phi_pi * kappa + phi_x + sigma - rho_u * (kappa + beta * phi_x + sigma * (1 + beta * (1 - rho_u)))));
checkdiff(g_yu, - (phi_pi - rho_u) / ((phi_pi - rho_u) * kappa + (phi_x + sigma * (1 - rho_u)) * (1 - beta * rho_u)));

if doTaylorWithNaturalRate
    checkdiff(g_piy, 0);
    checkdiff(g_yy, 1);
else
    checkdiff(g_piy, kappa / (1 - beta * rho_y) * (g_yy - 1));
    checkdiff(g_piy, -(kappa * sigma * (1 - rho_y)) / ((phi_pi - rho_y) * kappa + (phi_x + sigma * (1 - rho_y)) * (1 - beta * rho_y)) );
    checkdiff(g_yy, ((phi_pi - rho_y) * kappa + phi_x * (1 - beta * rho_y)) / ((phi_pi - rho_y) * kappa + (phi_x + sigma * (1 - rho_y)) * (1 - beta * rho_y)) );
end
%% effective Taylor system
Kx    = [g_yu (g_yy - 1)] / G;
Kybar = [0 1] / G;

if doTaylorWithNaturalRate
    barphi_pi = phi_pi + phi_x * Kx(1) + sigma * (rho_y - 1) * Kybar(1) ;
    barphi_y  = phi_x  * Kx(2) + sigma * (rho_y - 1) * Kybar(2) ;
else
    barphi_pi = phi_pi + phi_x * Kx(1) ;
    barphi_y  = phi_x  * Kx(2) ;
end

Abar = A;
Bbar = [1, -kappa; barphi_pi, barphi_y+sigma];
Dbar = [-1 kappa; 0 0];


BBbar      = Abar \ Bbar;
DDbar      = Abar \ Dbar;
lambdavals = eig(BBbar);


% check that original MSV still applies
checkdiff(Abar * G * F, Bbar * G + Dbar);
checkdiff(G * F, BBbar * G + DDbar);

GG    = reshape(-(kron(I,B) - kron(F',A)) \ D(:), Ny, Ny);
checkdiff(GG, G);


%% compute indeterminate outcomes
OmegaEpsilon = diag([sig_u, sig_y]);

[ZZ, TT] = schur(BBbar);
[ZZ, TT] = ordschur(ZZ,TT, 'udi');

Nchi       = 1;


ZZprime    = ZZ';
ZZprimeChi = ZZprime(:,1:Nchi);

TT11  = TT(1:Nchi, 1:Nchi) ;
invGZ  = GG \ ZZprimeChi;
FinvGZ = F * invGZ;
PSIgg = invGZ' / OmegaEpsilon * invGZ ;
PSIgf = invGZ' / OmegaEpsilon * FinvGZ ;
PSIff = FinvGZ' / OmegaEpsilon * FinvGZ;


b = 1 - TT11^2;
a = TT11^2 * PSIff - 2 * TT11 * PSIgf + PSIgg;
GammaBmax = sqrt(.25 * b^2 / a);


beliefscale_vec=0:0.0005:1;

for bb=1:length(beliefscale_vec)
beliefscale    = beliefscale_vec(bb);

GammaB      = beliefscale * GammaBmax;
SIGMAplus   = .5 * b / a + sqrt(.25 * b^2 / a^2  - GammaB^2 / a); 
SIGMAminus  = .5 * b / a - sqrt(.25 * b^2 / a^2  - GammaB^2 / a); 

minus_real=0;
plus_real=0;
if isreal(SIGMAplus)
SIGMA = SIGMAplus;
plus_real=1;
GammaVplus  = (- (SIGMA * ZZprimeChi') / (GG')  + TT11 * (SIGMA * ZZprimeChi') / (GG') * F') / OmegaEpsilon;
end

if isreal(SIGMAminus)
    minus_real=1;
SIGMA = SIGMAminus;
GammaVminus = (- (SIGMA * ZZprimeChi') / (GG')  + TT11 * (SIGMA * ZZprimeChi') / (GG') * F') / OmegaEpsilon;
end

%% out all together to create IRFs
% states are v(t) and delta(t)
% s(t+1) = AAA s(t) + BBB shocks(t+1)
% y(t)   = CCC s(t)

irflags = 16;

AAA = blkdiag(F, TT11);
if plus_real
BBBplus  = [eye(Ny, Ny + Nchi); GammaVplus GammaB] * blkdiag(chol(OmegaEpsilon, 'lower'), eye(Nchi));
end

if minus_real
BBBminus = [eye(Ny, Ny + Nchi); GammaVminus GammaB] * blkdiag(chol(OmegaEpsilon, 'lower'), eye(Nchi));
end
CCC = [GG ZZprimeChi];
CC  = [GG zeros(Ny,Nchi)];

if plus_real
[yirf, xirf]  = model2irf(AAA,BBBplus,CCC,irflags);
end

if minus_real
yirf_minus    = model2irf(AAA,BBBminus,CCC,irflags);
end
yirf_fullinfo = model2irf(AAA,BBBplus,CC,irflags);

ynames     = {'\pi', 'y'};
shocknames = {'cost push', 'natural rate', '\eta'};

% figure;
% plotirf(yirf, ynames, shocknames, yirf_fullinfo);
% sgtitle('Full vs Imperfect Info')
% 
% figure;
% plotirf(yirf, ynames, shocknames, yirf_minus);
% sgtitle('Imperfect Info: plus vs minus solution')
if bb==1
figure;
end
if plus_real
plotirf(yirf_fullinfo, ynames, shocknames,yirf);
end
hold on
if minus_real
plotirf(yirf_fullinfo, ynames, shocknames,yirf_minus);
end
hold on
sgtitle(sprintf('Full (blue) vs Imperfect Info '))

end

plotirf(yirf_fullinfo, ynames, shocknames);
hold on

print -depsc

savefig('figure_1.fig')
toc