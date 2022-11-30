function [Sigma, K, converged, iter, history] = abckalmani(A,B,C,K0,tolerance,maxiter)
% ABCKALMANI solves the Kalman filter using policy improvement iterations
% 
% The algorithm allows for singular in measurment innovations 
% (i.e. the rpocedure checks for the conditioning of VarZ and calls pinv
% when necessary)
%
% Usage: [Sigma, K, converged, iter, history] = abckalmani(A,B,C,K0,tolerance,maxiter)
%
% See also abckalman, abckalmanvfi, dare, pinv

%   Coded by  Elmar Mertens, em@elmarmertens.com

%% VERSION INFO
% AUTHOR    : Elmar Mertens
% $DATE     : 04-May-2010 11:45:27 $
% $Revision : 1.00 $
% DEVELOPED : 7.8.0.347 (R2009a)
% FILENAME  : abckalmani.m

 
if nargin < 4 || isempty(K0)
   K0 = zeros(size(C')); 
end

if nargin < 5 || isempty(tolerance)
   tolerance	= 1e-12;
end

if nargin < 6 || isempty(maxiter)
   maxiter     = 1e4;
end


if nargout > 4
   history = repmat(struct('K', NaN(size(C')), 'Sigma', NaN(size(A))), maxiter, 1);
else
   history = [];
end

BB          = B * B';
I           = eye(size(A));
iter        = 0;
converged   = false;

while ~converged && iter < maxiter
   iter = iter + 1;
     
   Sigma  = dlyapdoubling(A * (I - K0 * C), BB);
   
   covxz   = Sigma * C';
   varz    = C * Sigma * C';
   if rcond(varz) > 1e-10
      K = covxz / varz;
   else
      K = covxz * pinv(varz);
   end
   
   dK          = K - K0;
   converged   = max(abs(dK(:))) < tolerance;
   
   K0       = K;
   
   if ~isempty(history)
      history(iter).Sigma  = Sigma;
      history(iter).K      = K;
   end
   
end

if ~isempty(history)
   history = history(1:iter);
end

if ~converged
    if isriccatistable(A,B,C);
        error('em:msg', 'not converged, mae is %f\n\t But stability conditions have been met.', max(abs(dK(:))))
    else
        error('em:msg', 'not converged, mae is %f\n\t And stability conditions have NOT been met.', max(abs(dK(:))))
    end
end

Sigma = choppy(Sigma, -log10(tolerance));
K     = choppy(K, -log10(tolerance));


