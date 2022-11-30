function distance = integral_objective(D_vec,A,B,C,K,nobs,nshocks)

% Nx = size(A,1);
% Nw = size(B,2);

D=reshape(D_vec,nobs,nshocks);

BmKD     = B - K * D;
Omega    = BmKD * BmKD';

AmKC     = A - K * C;

% AmKCpowj = eye(size(A));
%
% ABBA = zeros(size(A));
%
% Sigma = zeros(Nx);

% version 1: call Matlab's dlyap
% Sigma    = dlyap(AmKC, Omega);
% version 2: call Hansen-Sargent doubling -- which is faster than dlyap!!
Sigma  = dlyapdoubling(AmKC, Omega);
% checkdiff(Sigma, Sigma2,1e-6);

FOC      = BmKD * D' + AmKC * Sigma * C';
distance = sum(FOC(:).^2);


end

