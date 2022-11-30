%% MODEL SETUP for Fisher-101

modellabel = 'fisher101';

%% calibrate parameters
horizon_simul=20000;
% real rate process
rho_r   = 0.9;
sigma_r = 1; % vol of r shocks



% policy rule parameters
phi_i  = 0;

phi_r  = 0;

%% setup model matrices
Nx   = 2;
Ny   = 1;
Nz   = 1;
Ni   = 1;
Neps = 2; % number of exogenous shocks
Nb   = Ny; % number of belief shocks
Xlabels   = {'\epsilon', '\nu^{\pi}'};
Ylabels   = {'\pi'};
YYlabels  = {'\pi', 'i'}; % YY is mathcal{Y}
EpsilonLabels = Xlabels(1:Neps);

ShockLabels = cat(2, EpsilonLabels, {'\eta^{\pi}'});


Ns = Nx + Ny;
Nw = Neps + Nb;

% setup matrices (as described in notes)
J = [1 0 0 
     0 1 0 
     0 0 1 ];

hatJ = zeros(Ns);

H = [rho_r 0 0 
     0   0 0 
     -1  0 0 ];

hatH = zeros(Ns);

Hi       = zeros(Ns,Ni);
Hi(Ns,1) = 1;

Bxeps  = diag([sigma_r, sigma_pi]);
BBx    = [Bxeps; zeros(Nx-Neps,Neps)];

% measurement
Cx = zeros(Nz,Nx);
Cx(1,2) = sigma_pi;

Cy = eye(Nz,Ny);   % C_z = I is a fixed assumption throughout
C  = [Cx Cy];

% policy rule
PHI_i = phi_i;
PHI_J = zeros(Ni,Ns);
PHI_H = [phi_r 0 phi_pi];
