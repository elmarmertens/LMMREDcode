%% Simple NK example from schwanLinREsetup.pdf

clear
clc
close all
tic
addpath ./../matlabbox/emtools/
addpath ./../schwantoolbox/

warning off

rng('default') % init random number generator

useDare = true;

%% Options for solution algorithm
tries_normal = 1000;
tries_uniform= 1000;
tries=tries_normal+tries_uniform; %number of candidate values
proposalLB=-500;
proposalUB=500; %upper and lower bound for uniform proposal
horizon_for_sol=1000; %truncation for computing infinite sum
fminTolFun = 1e-8; % tolerance for minimizer

fminOptions = optimset('Display', 'off','MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',fminTolFun,'TolX',1e-15);


%% other parameters

% IRF parameters
irflags  = 16;
aclags   = 20;
ticksize = 4;
%rng(0) % init random number generator

%% Load model setup


noisescale = 1;

setup301

%% full info solution

% JJ and AA are mathcal J and A in equation 20 of the notes:

JJ = [eye(Ni) zeros(Ni,Ns) zeros(Ni); ...
    zeros(Ns,Ni), J + hatJ, zeros(Ns,Ni); ...
    zeros(Ni) PHI_J zeros(Ni)];

AA = [zeros(Ni) zeros(Ni,Ns) eye(Ni); ...
    zeros(Ns,Ni) ,H + hatH, Hi; ...
    -PHI_i -PHI_H eye(Ni)];

Nxx  = Nx + Ni;
Nyy  = Ny + Ni;

% call Klein's solver -- linreABC is variant of Klein's solab
% maybe useful: linreABC uses matlab's own qz ordering rather than Sim's
% qzdiv
[GG, PP] = linreABC(JJ,AA,Nxx);

% same results with Klein's own routine:
% [checkGG, checkPP] = solab(JJ,AA,Nxx);
% checkdiff(checkGG, GG);
% checkdiff(checkPP, PP);



BB       = [zeros(Ni,Neps); BBx];

PPi    = PP(1:Ni,:);
GGy    = GG(1:Ny,:);
II     = [zeros(Nx,Ni) eye(Nx)];
FF     = [II; GGy];
GGyx   = GGy(:,Ni+(1:Nx));
GGyeps = GGy * BB; % same as GGyx * BBx;
checkdiff(GGyeps, GGyx * BBx);

% plot some irf
% yirf_fullinfo = zeros(Nyy, Nw, irflags+1);
% yirf_fullinfo(:,1:Neps,:) = model2irf(PP,BB, GG, irflags);
% plotirf(yirf_fullinfo(:,1:Neps,:), YYlabels, EpsilonLabels)

%% check full info solution
switch lower(modellabel)
    case 'fisher101'
        gbar = GGyx(2);
        checkdiff(gbar, 1 / (phi_pi - rho_r));
end

%% set up residual system
A      = J \ H;
hatAsx = J \ (hatH * FF - hatJ * FF * PP + Hi * PPi);

Axx = A(1:Nx,1:Nx);
Axy = A(1:Nx,Nx+(1:Ny));
Ayx = A(Nx+(1:Ny),1:Nx);
Ayy = A(Nx+(1:Ny),Nx+(1:Ny));
Gx  = eye(Nx,Ns);

LLL = Cx+GGyx;

%setting up matrices for ABCD system


AAA = Axx-Axy*Cx;

BBB = [BBx zeros(size(BBx,1),Nb)];

CCC = Cx*(Axx-Axy*Cx)+Ayx-Ayy*Cx;

%% setup ABCD and call numerical solver


% allocate memory for results from tries
solutions=NaN(tries,(Neps+Nb)*Ny);
objective_vec=NaN(tries,1);
K_array=NaN(Nx,Ny,tries);
sol_ind = NaN(tries,1);

% allocate storage for initial value
D0_vec = NaN(tries,Nz * Nw);



% loop over tries


for tt=1:tries
  
    % setup solver
    obj = @(D_vec) (riccati_objective(D_vec,AAA,BBB,CCC,LLL,Ny,Nw,useDare));
    
    if tt <= tries_normal
        D0_vec(tt,:) = randn(1,Nz * Nw);
    else
        D0_vec(tt,:) = 5 + 10 * rand(1,Nz * Nw);
    end

    % call solver

    try
        [solutions(tt,:),objective_vec(tt),exitflag, fminOutput] = fminunc(obj, D0_vec(tt,:), fminOptions);
    catch
        exitflag =-10;
        objective_vec(tt)=10;

    end
    if exitflag < 1 || objective_vec(tt) > fminTolFun
        sol_ind(tt)=0;
    else
        sol_ind(tt)=1;

        % collect Kalman gain
        [~,KKK] = obj(solutions(tt,:));

        K_array(:,:,tt)=KKK;
       
    end
    
    % progressbar(tt/tries)
end


nsol = sum(sol_ind);
DDD  = reshape(solutions(logical(sol_ind),:)', Ny, (Neps+Nb), nsol);

fprintf('Found %d solutions after %d tries (%d with normal proposals, %d uniform).\n', nsol, tries, tries_normal, tries_uniform)
fprintf('Lowest value for objective was %6.4e.\n', min(objective_vec))



%% plot results

% generate full info IRF
yirf_fullinfo = zeros(Nyy, Nw, irflags+1);
yirf_fullinfo(:,1:Neps,:) = model2irf(PP,BB, GG, irflags);

yirf_simul_full_info=modelsimul(PP,BB, GG, horizon_simul);

% collect partial info
K_array_final=K_array(:,:,logical(sol_ind));
if isempty(K_array_final)
    error('no joy: could not find any solution!')
end


yirf_simul       = NaN(Nyy, horizon_simul + 1, nsol); 
yirf_partialinfo = NaN(Nyy, Nw, irflags + 1, nsol);
    
parfor ind_loop=1 : size(K_array_final,3)
    
    [ A_final, B_final,C_final,Chat_final, H_final, Hhat_final] = state_space_final( A,BBx,C,DDD(:,:,ind_loop),K_array_final(:,:,ind_loop),PP,GG,GGyx,Nx,Ny );
    
    
    yirf_partialinfo(:,:,:,ind_loop)  = model2irf(A_final,B_final, C_final, irflags);
    yirf_simul(:,:,ind_loop)          = modelsimul(A_final,B_final, C_final, horizon_simul);
     
    
end

% plot IRF
showlags = 16;
figure
hold on
for s = 1 : nsol
    plotirf(yirf_fullinfo(:,:,1:showlags+1), YYlabels, ShockLabels,yirf_partialinfo(:,:,1:showlags+1, s), 4)
end
plotirf(yirf_fullinfo(:,:,1:showlags+1), YYlabels, ShockLabels);

print('figure_8', '-depsc')
savefig('figure_8.fig')

plot_simul_results2em;

toc

