function [distance,K] = riccati_objective(D_vec,A,B,C,L,nobs,nshocks, useDare)

if nargin < 7
    useDare = true; % much quicker
end

% Nx = size(A,1);
% Nw = size(B,2);

% solve riccati equation to obatin Kalman gain given D
% then check whether L * K = I

D=reshape(D_vec,nobs,nshocks);


% [Sigma, ~, K] = dare(A', C', B*B', D * D', B * D');
% K = K'; % adapting from dare notation
% checkdiff(K, (A * Sigma * C' + B * D') / (C * Sigma * C' + D * D'));

% via Matlab's dare (solving matrix pencil)
if useDare
    try % dare may break if parameters violate detectability etc.
        [~, ~, K] = dare(A', C', B*B', D * D', B * D');
        K = K'; % adapting from dare notation
    catch
        distance = 1e2; % some penalty
        return
    end
else
    % via VFI riccati
    [~, F] = riccati(A', C', B*B', D * D', B * D');
    % checkdiff(-F',K);
    K = -F'; % adapting from riccati notation
end

FOC      = L * K - eye(nobs);
distance = sum(FOC(:).^2);


end

