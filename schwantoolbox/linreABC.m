function [f,p,n,l,M] = linreABC(a, b, nk , c, PHI, PSI)
% function [f,p,n,l,M,abcd] = linreABC(a, b, nk , c, PHI, PSI)
% Purpose: Solves for the recursive representation of the stable solution to a system
% of linear difference equations. And allows for exogenous drivers with VAR dynamics (companion matrix PHI)
%
% Inputs: Two square matrices a and b, natural number nk, matrices C and PHI (square)
%
% a, b and c are the coefficient matrices of the difference equation
%
% a*x(t+1|t) = b*x(t) + c * z(t)
% z(t+1)     = PHI * z(t) + PSI w(t) 
%
% note: PSI is irrelevant for the solution, matters only for putting together the abcd state space
% 
% where x(t) is arranged so that the state variables come first, and
% nk is the number of state variables.
%
% Outputs: the decision rule f and the law of motion p. If we write
%
% x(t) = [k(t);u(t)] where k(t) contains precisely the state variables, then
% 
% u(t)     = f*k(t) + n * z(t) 
%
% and
%
% k(t+1|t) = p*k(t) + l * z(t)
%
% See also: solab, solab2, solabc, solaborg, qzdiv, qzswitch
%
% adapted from Paul Klein's solab.m by Elmar Mertens
% this version replaces Sims' qzdiv with Matlab's ordqz
% www.elmarmertens.ch

%   Coded by  Elmar Mertens, em@elmarmertens.com

nu  = size(a,1) - nk;

if nargin < 4
   c     = [];
   PHI   = [];
   PSI   = [];
end
nz    = size(c, 2);
if nargin < 6
   PSI   = eye(nz);
end
nw  = size(PSI, 2);

[s,t,q,z] = qz(a,b);  % factorization of the matrix pencil b-za, see below for ordering


% em: generalized eigenvalues and check for regularity of pencil
sii            = diag(s);
tii            = diag(t);
nulli          = abs(sii) < eps;

% compute generalized eigenvalues
geigenvalues            = Inf(size(a,1), 1); 
geigenvalues(~nulli)    = abs(tii(~nulli)) ./ abs(sii(~nulli));

orgeig         = geigenvalues;
geigenvalues           = choppy(geigenvalues, 12);
if ~isequal(geigenvalues <= 1, orgeig <= 1)
   warning('em:msg', 'there were some eigenvalues just on the boundary of the unit circle, chopping off at 1e-12')
end

Nstable  = sum(geigenvalues <= 1);
if any(geigenvalues == 1)
   warning('em:unitroot', 'the system has %d unit root(s) (will be treated as "stable")', sum(geigenvalues == 1))
end

if ~isempty(nk) && Nstable ~= nk
   error('em:msg', 'there are %d stable eigenvalues but you set nk=%d', Nstable, nk)
end

% re-order QZ decomposition
unstable        = geigenvalues > 1;
[s,t,q,z]       = ordqz(s,t,q,z, ~unstable);
sii             = diag(s);
tii             = diag(t);

if any(~any(abs([sii tii]) > 1e-10, 2)) 
   warning('em:msg', 'not a regular pencil')
end

% % correct number of stable eigenvectors
% if abs(t(nk,nk))>abs(s(nk,nk)) || abs(t(nk+1,nk+1))<abs(s(nk+1,nk+1));
%    warning('Wrong number of stable eigenvalues.');
% end

z11 = z(1:Nstable,1:Nstable);
z12 = z(1:Nstable,Nstable+1:end);
z21 = z(Nstable+1:end,1:Nstable);
z22 = z(1+Nstable:end,1+Nstable:end);
s22 = s(1+Nstable:end,1+Nstable:end);
s12 = s(1:Nstable,Nstable+1:end);

t12 = t(1:Nstable,Nstable+1:end);
t22 = t(1+Nstable:end,1+Nstable:end);

% z12 = z(1:Nstable,Nstable+1:end);
% z22 = z(1+Nstable:end,1+Nstable:end);
% s22 = s(1+Nstable:end,1+Nstable:end);
% s12 = s(1:Nstable,Nstable+1:end);
% t12 = t(1:Nstable,Nstable+1:end);
% t22 = t(1+Nstable:end,1+Nstable:end);
% nu  = size(a,1) - Nstable;

if rank(z11) < Nstable
   error('em:msg', 'Invertibility condition violated: rank(z11) = %d < Nstable = %d.', rank(z11), Nstable)
   %    fprintf('By the way, since you are running into this non-invertbility error: sort(geigenvalues)\n')
   %    fprintf('%20.4e \n', sort(geigenvalues))
   %    [f,p,n,l] = deal([]);
   %    status.code = false;
   %    status.msg = sprintf('Invertibility condition violated: rank(z11) = %d < Nstable = %d.', rank(z11), Nstable);
   %    return
end

z11i  = z11 \ eye(Nstable);
s11   = s(1:Nstable,1:Nstable);
t11   = t(1:Nstable,1:Nstable);

f = z21 * z11i;
p = z11 / s11 * t11 * z11i;

%% contruct responses to exogenous states
if ~isempty(c)
   
   %% construct M
   q2c   = q(nk+1:end,:) * c;
   Inz   = eye(nz);
   
   % recursion
   M     = zeros(nu, nz);
   for i = nu : -1 : 1
      rip = q2c(i,:);
      for j = i + 1 : nu
         rip = rip + (t22(i,j) * M(j,:) - s22(i,j) * M(j,:) * PHI);
      end
      M(i,:) = rip / (s22(i,i) * PHI - t22(i,i) * Inz);
   end
   
   % some slower alternatives (slower even for moderate systems)
   % M     = reshape((kron(PHI', s22) - kron(Inz, t22)) \ q2c(:), nu, nz);
   % M = dlyap(t22 \ s22,  PHI, -t22 \ q2c)
   
   
   %% L, N
   z12M  = z12 * M;
   n     = z22 * M - f * z12M;
   l     = - p * z12M + z11 / s11 * (t12 * M - s12 * M * PHI + q(1:nk, :) * c) + z12M * PHI;
   
else
   [n, l, M] = deal([]);
end

%% take out imaginary parts

if any(max(abs(imag([f(:); p(:); n(:); l(:)]))) > 1e-10)
   warning('em:msg', 'Some imaginary parts are left (even though they shouldn''t)\n\tlargest element has size %e', max(abs(imag([f(:); p(:); n(:); l(:)]))))
end

f = real(f);
p = real(p);
n = real(n);
l = real(l);

%% transform into state space
% if nargout > 5
%    abcd.A   = [PHI zeros(nz, nk); l p];
%    abcd.B   = [PSI; zeros(nk, nw)];
%    abcd.Cx  = [n f];
%    abcd.C   = abcd.Cx * abcd.A;
%    abcd.D   = abcd.Cx * abcd.B;
%    abcd = abcddims(abcd);
% end

%% check solution
% policymatrix = [eye(Nstable); f];
% if checkdiff(a * policymatrix * p, b * policymatrix, [], 'LinDiff Check')
%     warning('check lindiff')
%     keyboard
% end


function x = choppy(x, d)
% function x = choppy(x, d)
% rounds up to d'th digit after decimal
% difference with chop: chop retains d first *significant* (non-zero) digits
% default: d = 12
%
% Alternartive usages: 1) d is positive integer
%                      2) when d < 1 (but positive), it will be interpreted
%                      like a tolerance, say 1e-10, i.e. x will be cut to
%                      the same amount of digits as d

%   Coded by  Elmar Mertens, em@elmarmertens.com

if nargin < 2
   d = 12;
end

if d >= 1
   d = ceil(d);
   x = round(x * (10^d)) / 10^d;
elseif d > 0
   x = round(x / d) * d;
else
   error('argument d must be positive (here: d=%e)', d)
end