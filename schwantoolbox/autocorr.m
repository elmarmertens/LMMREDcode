function acf = autocorr(y,numLags)
%AUTOCORR Sample autocorrelation
%
% Syntax:
%
%
%   [acf,lags] = autocorr(y,numLags)
%
% Poor man"s version of Matlab's autocorr


[rows,columns] = size(y);

if (rows ~= 1) && (columns ~= 1)
    
    error(message('econ:autocorr:NonVectorInput'))
      
end

rowSeries = (size(y,1) == 1);

y = y(:);         % Ensure a column vector
N = length(y);    % Sample size
defaultLags = 20; % Recommendation of [1]

% Ensure numLags is a positive integer or set default:

if (nargin >= 2) && ~isempty(numLags)
    
   if numel(numLags) > 1
       
      error(message('econ:autocorr:NonScalarLags'))
        
   end
   
   if (round(numLags) ~= numLags) || (numLags <= 0)
       
      error(message('econ:autocorr:NonPositiveInteger'))
        
   end
   
   if numLags > (N-1)
       
      error(message('econ:autocorr:LagsTooLarge'))
        
   end
   
else
    
   numLags = min(defaultLags,N-1); % Default
   
end



% Convolution, polynomial multiplication, and FIR digital filtering are all
% the same operation. The FILTER command could be used to compute the ACF
% (by convolving the de-meaned y with a flipped version of itself), but
% FFT-based computation is significantly faster for large data sets.

% The ACF computation is based on [1], pages 30-34, 188:

nFFT = 2^(nextpow2(length(y))+1);
F = fft(y-mean(y),nFFT);
F = F.*conj(F);
acf = ifft(F);
acf = acf(1:(numLags+1)); % Retain non-negative lags
acf = acf./acf(1); % Normalize
acf = real(acf);

