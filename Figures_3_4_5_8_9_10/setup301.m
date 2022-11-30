%% MODEL SETUP for NK201

modellabel = 'NK301';

% NK model as in 201; but tracking also the level of actual GDP

horizon_simul=20000;
%% calibrate parameters

% preferences
beta    = .99;
phi     = 1;   % inverse Frisch labor elasticity
sigma   = 1;   % Risk-Aversion / Inverse IES

% Phillips Curve
gamma   = 0.25;   % weight on bwd-looking inflation in PC
theta   = .75;    % Calvo Probability of *not* resetting price
kappa   = (1 - theta) * (1 - beta * theta) / theta * (sigma + phi); 

% real rate process
rho_dybar   = 0.75;
sigma_dybar = .001*400*.75; % vol of r shocks %values roughly in line with our JME
% we know that r_t=sigma E_t delta(ybar_t+1)=sigma*rho*delta(y_bar_t)
%Suppose we have data on the measurement error in GDP:
%delta(y_obs)=delta(y)+meas where meas could be persistent
%If this is consistent with a mismeasured gdp level y_obs=y+meas_level
%(Where meas_level is iid) the moments of meas should be equal to those of
%2*meas_level. In the JME with Thomas we estimated meas to have
%persistence of 0.5 and a standard deviation of .006. Inflation has
%persistence of .1 and standard deviation of .002.

pi_vol_uncond=sqrt(.002^2/(1-.1^2))*400;
y_vol_uncond=sqrt(.006^2/(1-.5^2))*400*.5;

% scale factor to scale noise up/down
if ~exist('noisescale', 'var')
    noisescale = 1;
end

% noise volatiltiies
sigma_pi = pi_vol_uncond * noisescale;
sigma_y  = y_vol_uncond * noisescale;

% policy rule parameters
phi_i  = 0;
phi_pi = 2.5;
phi_x  = 0.5;
phi_r  = 0;

%% setup model matrices
Nx   = 5;
Ny   = 2;
Nz   = 2;
Ni   = 1;
Neps = 3; % number of exogenous shocks
Nb   = Ny; % number of belief shocks
Xlabels   = {'\Delta y', '\nu^{\pi}', '\nu^x', 'pi(-1)', 'ybar(-1)'};
Ylabels   = {'pi', 'x'};
YYlabels  = {'\pi', 'x', 'i'}; % YY is mathcal{Y}
EpsilonLabels = Xlabels(1:Neps);

ShockLabels = cat(2, EpsilonLabels, {'\eta^{\pi}'},{'\eta^x'});


Ns = Nx + Ny;
Nw = Neps + Nb;

% setup matrices (as described in notes)
J = [1 0 0 0   0              0    0
     0 1 0 0   0              0    0
     0 0 1 0   0              0    0
     0 0 0 1   0              0    0
     0 0 0 0 1 0 0
     0 0 0 (1- gamma *beta) 0 -beta 0
     0 0 0 0                0 1     sigma];

hatJ = zeros(Ns);

H = [rho_dybar 0 0 0 0     0    0
     0   0 0 0 0    0    0
     0   0 0 0 0    0    0
     0   0 0 0 0    1    0
     1 0 0 0 1 0 0
     0   0 0 gamma 0 0 kappa 
     -sigma * rho_dybar  0 0 0     0 0 sigma];

hatH = zeros(Ns);

Hi       = zeros(Ns,Ni);
Hi(Ns,1) = 1;

Bxeps  = diag([sigma * rho_dybar * sigma_dybar, sigma_pi, sigma_y]);
BBx    = [Bxeps; zeros(Nx-Neps,Neps)];

% measurement
Cx = zeros(Nz,Nx);
Cx(1,2) = sigma_pi;
Cx(2,3) = sigma_y;
Cx(2,1) = 1;
Cx(2,5) = 1;


Cy = eye(Nz,Ny);   % C_z = I is a fixed assumption throughout
C  = [Cx Cy];

% policy rule
PHI_i = phi_i;
PHI_J = zeros(Ni,Ns);
PHI_H = [phi_r * sigma * rho_dybar 0 0 0 0 phi_pi phi_x];
