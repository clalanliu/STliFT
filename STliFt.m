%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Semidefinite relaxation-based algorithm      %
%               with MATLAB Implementation             %
%                                                      %
% Author: Chang Le Liu                      07/06/2018 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This code is the Matlab implimentation of K. Jaganathan et al 
% "STFT Phase Retrieval: Uniqueness Guarantees and Recovery Algorithms,"
% IEEE Journal of Selected Topics in Signal Processing ( Volume: 10, Issue: 4, June 2016 )

function x_restore = STliFt(S, xlen, winLen, hop, N, sflag)
% function: x_restore = STliFt(S, xlen, N, R)
% S - STFT magnitude measurements
% xlen - Length of the original signal producing S
% winLen - Size of the hamming symmetric window
% hop - Seperation in time between adjacent short-time sectoins

% Optinal args
% N - Number of DFT points. Default value is size(S,1) 
% R - Number of short-time section considered. Default value is size(S,2) 
% sflag - Window sampling, specified as one of the following:
%   'symmetric'(default) - X Use this option when using windows for filter design.
%                        - (Used in spectrogram() in MATLAB)
%   'periodic'           - X This option is useful for spectral analysis because it
%                          enables a windowed signal to have the perfect periodic 
%                          extension implicit in the discrete Fourier transform. 
%                          When 'periodic' is specified, hamming computes a window
%                          of length L + 1 and returns the first L points.


%% parameters
if nargin < 5
    N = size(S,1);
end
if nargin < 6
    sflag = 'symmetric';
end

S = S.^2;

% magnitude-square
R = size(S,2);

%% N-point DFT matrices
D = dftmtx(N);
Ds = cell(N,1);
parfor m=1:N
    Ds{m,1} = circshift(D(m,:),(m-1)*hop)'*circshift(D(m,:),(m-1)*hop);
end

%% W matrix array
win = hamming(winLen, sflag);
winPad = padarray(win, N-winLen,'post');
Ws = cell(R,1);
parfor r=1:R
    Ws{r,1} = diag(circshift(winPad,(r-1)*hop));
end

%% CVX SDP optimazation
cvx_begin sdp
    variable X(N,N) semidefinite
    minimize(trace(X))
    for m=1:N
        for r=1:R
            trace(Ws{r,1}'*Ds{m,1}*Ws{r,1}*X) == S(m,r)
        end
    end
cvx_end

%%
X_full = full(X);
[~, x_restore] = bestRankOneApproximation(X_full);

end

function [apprA, apprU] = bestRankOneApproximation(A)
    [U,S,V]=svd(A);
    u = U(:,1);
    v = V(:,1);
    apprU = -sqrt(S(1,1))*u;
    apprA = S(1,1)*u*v';
end