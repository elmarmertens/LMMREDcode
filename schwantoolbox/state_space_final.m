function [ A_final, B_final,C_final, Chat_final, H_final, Hhat_final] = state_space_final( A,BBx,C,DDD,K_x,PP,GG,GGyx,Nx,Ny )
%Uthis code sets up the state space system for the actual variables. The
%observable vector is ordered as [X_t' Y_t' i_t]'
%The state equation is given by equation (50) in the linREsetup notes,
%whereas the measurement equation is the combination of the expressions
%just below equation (50)

% EM: amended from earlier version by CM Dec 18 2017:
% - dropped D_final from list of outputs (always zero anyway)
% - Y vector created by C_final contains Y and i (as in fullinfo version of the model)
% - Chat_final creates central bank projections of Y
% - H_final and Hhat_final generate X and projections of X

Ni = 1; % currently assumed

Nz  = Ny;
Ns  = Nx + Ny;
Nxx = Nx + Ni;
Nw  = size(DDD,2);
Neps = size(BBx,2);
Nb  = Nw - Neps;

K_y=GGyx*K_x;
K=[K_x;K_y];
KK_x=[zeros(Ni,Nz);K_x];

Cx      = C(:,1:Nx);
B_etab   = DDD(:,Neps+1:end);
B_etaeps = DDD(:,1:Neps) - Cx * BBx; % need to strip out shocks to Z induced by Cx * X*(t)

B=[BBx zeros(Nx,Nb); B_etaeps B_etab];

A_final=[A*(eye(Ns)-K*C) zeros(Ns,Nxx);PP*KK_x*C PP];
B_final=[B;zeros(size(PP,1),Nw)];

% C_final=[eye(Nx) zeros(Nx,Ny+1) eye(Nx);zeros(Ny,Nx) eye(Ny) GG(1:Ny,:);PP(1,:)*KK_x*C  PP(1,:)];
C_final    = [zeros(Ny,Nx) eye(Ny) GG(1:Ny,:);PP(1,:)*KK_x*C  PP(1,:)];
Chat_final = [K_y * C GG(1:Ny,:);PP(1,:)*KK_x*C  PP(1,:)];

H_final    = [eye(Nx) zeros(Nx,Ny+1) eye(Nx)];
Hhat_final = [K_x * C zeros(Nx,1) eye(Nx)];


end
