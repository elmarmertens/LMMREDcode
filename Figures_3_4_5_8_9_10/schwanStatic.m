function B_eta = schwanStatic(Ceps, GGyeps, solutionType)
% SCHWANSTATIC ... 
% simple solution, mark 1 
%   ... 

%% VERSION INFO 
% AUTHOR    : Elmar Mertens 
% $DATE     : 05-Oct-2017 13:49:05 $ 
% $Revision : 1.00 $ 
% DEVELOPED : 9.3.0.713579 (R2017b) 
% FILENAME  : schwanStatic.m 

if nargin < 3
    solutionType = 1;
end

[Nz, Neps]  = size(Ceps); %#ok<ASGLU>
[Ny, check] = size(GGyeps);

if check ~= Neps
    error('dimension mismatch: Neps=%d, check=%d', Neps, check)
end
Ne = Neps + Ny; %#ok<NASGU>



% map into call to staticInfoPolicy
% solutionType = 0;
alpha = 1;
beta  = 1;

[Betaeps, Betab] = staticInfoPolicy(GGyeps, eye(Neps), Ceps, solutionType, [], alpha, beta);

B_eta   = [Betaeps Betab];
B_eta   = choppy(B_eta); % rounding things to 12th digit