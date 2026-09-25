function CRPS=CRPSg(xmean,xstd,yref)
% CRPS=CRPSg(xmean,xstd,yref)
% A Gaussian approximation of the Continuous Ranked Probability Score
% Input
%   xmean : Mean of the estimated (typically posterior) distribution
%   xstd  : Std of the estimated distribution
%   yerf  : Reference turth value that we are comparing to
% Output
%   CRPS  : The estimated CRPS
% All inputs and outputs are all typically either scalars or column vectors
% with corresponding entries (e.g. different experiments or variables)
% Based on the Gaussian approximation from Gneiting https://doi.org/10.1175/MWR2904.1
deg=xstd<eps(1); % Check if the distribution is degenerate
err=yref-xmean; % Error-like term
z=err./xstd; % z-score
phi=normpdf(z);
Phi=normcdf(z);
CRPS=zeros(size(err));
CRPS(deg)=abs(err(deg)); % Reduces to absolute error in degenerate case
CRPS(~deg)=xstd(~deg).*(z(~deg).*(2.*Phi(~deg)-1)+2.*phi(~deg)-1/sqrt(pi));
% The above is eq. 5 in the paper of Gneiting
end
