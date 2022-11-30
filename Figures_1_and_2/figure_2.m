%% schwanNK example with inflation observed (and cost push shock only)

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

phi_pi = 1.5; 
phi_x  = 1;   %  .5 leads to barphi_pi = 1 (with phi_pi = 1.5)

% stochastic process for cost-push shocks
sig_u  = 1; % variance
rho_u  = 0;


Ny      = 2;
Nv      = 1;

%% pi-y full info system
% A y(t+1) = B y(t) + D v(t)
% y(t) = [pi, y]

A = [beta 0; 1 sigma];
B = [1, -kappa; phi_pi, phi_x+sigma];
D = [-1; 0];

F = rho_u;
I = eye(Nv);

Gvec = -(kron(I,B) - kron(F',A)) \ D(:);
G    = reshape(Gvec, Ny, Nv);


checkdiff(A * G * F, B * G + D);

g_piu = G(1,1);
g_yu  = G(2,1);

%% check coefficients
checkdiff(g_piu, (1 + kappa * g_yu) / (1 - beta * rho_u));
checkdiff(g_yu, - (phi_pi - rho_u) / (phi_pi * kappa + phi_x + sigma - rho_u * (kappa + beta * phi_x + sigma * (1 + beta * (1 - rho_u)))));
checkdiff(g_yu, - (phi_pi - rho_u) / ((phi_pi - rho_u) * kappa + (phi_x + sigma * (1 - rho_u)) * (1 - beta * rho_u)));

%% effective Taylor system
Kx    = g_yu / g_piu;

barphi_pi = phi_pi + phi_x * Kx(1);
checkdiff(barphi_pi, phi_pi * sigma / (phi_x + sigma));

Abar = A;
Bbar = [1, -kappa; barphi_pi, sigma];
Dbar = [-1; 0];


BBbar      = Abar \ Bbar;
DDbar      = Abar \ Dbar;
lambdavals = eig(BBbar);



% check that original MSV still applies
checkdiff(Abar * G * F, Bbar * G + Dbar);
checkdiff(G * F, BBbar * G + DDbar);

GG    = reshape(-(kron(I,B) - kron(F',A)) \ D(:), Ny, Nv);
checkdiff(GG, G);

%% compute indeterminate outcomes

ndxpi          = 1; % index into y to pin-point inflation
OmegaEpsilon   = sig_u;

[ZZ, TT] = schur(BBbar);
[ZZ, TT] = ordschur(ZZ,TT, 'udi');

Nchi       = 1;


ZZprime    = ZZ';
ZZprimeChi = ZZprime(:,1:Nchi);

TT11  = TT(1:Nchi, 1:Nchi) ;
invGZ  = GG(ndxpi,:) \ ZZprimeChi(ndxpi,:);
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

if isreal(SIGMAplus)
SIGMA = SIGMAplus;
plus_real=1;
GammaVplus  = (- (SIGMA * ZZprimeChi(ndxpi,:)') / (GG(ndxpi,:)')  + TT11 * (SIGMA * ZZprimeChi(ndxpi,:)') / (GG(ndxpi,:)') * F') / OmegaEpsilon;
end

if isreal(SIGMAminus)
    minus_real=1;
SIGMA = SIGMAminus;
GammaVminus = (- (SIGMA * ZZprimeChi(ndxpi,:)') / (GG(ndxpi,:)')  + TT11 * (SIGMA * ZZprimeChi(ndxpi,:)') / (GG(ndxpi,:)') * F') / OmegaEpsilon;
end

%% out all together to create IRFs
% states are v(t) and delta(t)
% s(t+1) = AAA s(t) + BBB shocks(t+1)
% y(t)   = CCC s(t)

irflags = 16;

AAA = blkdiag(F, TT11);

if plus_real
BBBplus  = [eye(Nv, Nv + Nchi); GammaVplus GammaB] * blkdiag(chol(OmegaEpsilon, 'lower'), eye(Nchi));
end

if minus_real
BBBminus = [eye(Nv, Nv + Nchi); GammaVminus GammaB] * blkdiag(chol(OmegaEpsilon, 'lower'), eye(Nchi));
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

ynames     = {'inflation', 'output'};
shocknames = {'cost push', 'belief'};

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

savefig('figure_2.fig')

toc

