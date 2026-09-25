function [w,logZ] = PIES( obs, pred, R, priormean, priorcov, proposal, propmean, propcov )
%% PIES w = PIES( obs, pred, R, priormean, priorcov, proposal )
%
% Weighting step for a Particle-adjusted Iterative Ensemble Smoother (PIES)
% Inputs:
%    obs: Observation vector (m x 1 array)
%    pred: Predicted observation ensemble matrix (m x N array)
%    R: Observation error covariance 'matrix' (m x 1 array, or scalar)
%    priormean: The mean (vector) of the prior (n x 1 array)
%    prirocov: The covariance (matrix) of the prior (n x n array)
%    proposal: Samples from the proposal (n x N array)
%        
%    Outputs:
%        w: Posterior weights (N x 1 array)
%    Dimensions
%        N is the number of ensemble members, n is the number of state variables and/or
%        parameters, and m is the number of observations.
%        
% The PIES scheme is obtained by using the output of the Iterative Ensemble Smoother, 
% i.e. a Gaussian distribution, as the proposal in importance sampling. Note that we do 
% not need to use the final posterior from IES as the proposal, we can instead use the
% "posterior" at any iteration l=0,...,Na where Na is the number of assimilation cycles.
% Note that the in the special case l=0 we end up with the standard PBS scheme if we disregard
% the differences between the true and sampled prior statistics.
%   
% Code by: K. Aalstad as part of the Spot-On project 
% The scheme is derived in https://doi.org/10.5194/amt-15-7293-2022

No=numel(obs);
Ne=size(pred,2);

if size(R,1)==1&&size(R,2)==1
    R=R.*eye(No);
elseif size(R,1)==1||size(R,2)==1
    R=diag(R);
end

A0=(proposal-priormean);
AT=A0';
B=priorcov\A0;
Phi0=-0.5.*sum(AT.*(B'),2); % Equivalent but faster than below for large Ne
propm=propmean; % Ns x Ne
A=proposal-propm; % Ns x Ne
C=propcov; % Ns x Ns
Cinv=pinv(C);
AT=A';
B=Cinv*A;
Phip=-0.5.*sum(AT.*(B'),2); % Equivalent but faster than below for large Ne
residual=obs-pred; % No x Ne
AT=residual';
B=R\residual;
Phid=-0.5.*sum(AT.*(B'),2);


Phi=Phid+Phi0-Phip;
Phimax=max(Phi);
Phis=Phi-Phimax;

w=exp(Phis);
logZ=logsumexp(Phi(:)); % Log evidence approximation (ignoring constants)
% Important not to include Phimax here since it's not constant!
w=w./sum(w);
if any(isnan(w))
    disp('nan weights in PIES')
end
end
