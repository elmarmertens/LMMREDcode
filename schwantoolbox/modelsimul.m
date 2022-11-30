function [yirf, xirf] = modelsimul(A,B,C,lags)
% MODEL2IRF
%  computes impulse responses for state space model with matrices A,B,C
%  simulates IRF for periods 0:lags with shock occuring at t=shocktime
%  outputs have dimension Ny x Nw x (lags+1) respectively Nx x Nw x (lags+1)
%
% USAGE: [yirf, xirf] = model2irf(A,B,C,lags,shocktime)
% defaults: lags = 16; shocktime = 0;

%   Coded by  Elmar Mertens, em@elmarmertens.com


error(abcdchk(A,B,C)) % check dimensional consistency

nx = size(A, 1);
nw = size(B, 2); 
ny = size(C, 1);


% init irf as 3D matrix
yirf = NaN(ny, lags + 1); 
xirf = NaN(nx, lags + 1); 
% notice: NaN is useful in partiuclar hen shocktime >0
% using zeros instead, might distort IRF plots (things will appear to "rise" before the shock, when connecting the dots
% (thanks to Bob King for pointing this out)

% draw random variables
shocks = randn(nw,lags + 1);

% impact responses
k = 0 ;  % notice: matlab starts indexing only at 1, our time index can start at 0
xirf(:,k + 1) = B * shocks(:,k + 1); % unit s.d. shocks
yirf(:,k + 1) = C * xirf(:,k + 1);

% simulate rest
for k = 1 : lags
   xirf(:,k + 1) = A * xirf(:,k) + B * shocks(:,k + 1);     % notice: A^(k-shocktime) * xirf(:,:,shocktime + 1) works too but is slower to compute
   yirf(:,k + 1) = C * xirf(:,k + 1);
end
