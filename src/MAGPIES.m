function [w,logZ] = MAGPIES( obs, pred, Rvec, prim, priv, pros, prom, proc )
%% MAGPIES w = MAGPIES( obs, pred, Rvec, prim, priv, pros, prom, proc )
%
% NB: This function is the same as AMIS with a Gaussian DMP obtained from
% several generations of the IES scheme (externally to this function).
%
% Weighting step for the Multiple Adaptive Guided Particle Iterative
% Ensemble Smoothers (MAGPIES) scheme.
%
% Dimensions: No = Number of observations in the batch to assimilate.
%             Np = Number of parameters to update.
%             Ne = Number of ensemble members (particles) per iteration
%             Na = Number of assimilation cycles in IES
%             Ng = Number of generations of IES
%             Nq = Number of Gaussian proposals = (Na+1) x Ng
%             Nh = Number of historical particles = Ne x Nq
% 
% Inputs:
%    obs: Observation vector (No x 1 array)
%    pred: Predicted observation ensemble matrix (No x Ne x Nq array)
%    Rvec: Observation error variance (No x 1 array)
%    prim: The mean vector of the prior (Np x 1 array)
%    priv: The variance vector of the prior (Np x 1 array)
%    pros: Samples from the proposal (Np x Ne x Nq arrray)
%    prom: Mean of the proposal (Np x Nq array)
%    proc: Covariance of the proposal (Np x Np x Nq array)
%        
%    Outputs:
%        w: Posterior weights (Nh x 1 array)
%        logZ : Log evidence (marginal likelihood)
%    Dimensions

%    The MAGPIES scheme is obtained by combining the PIES scheme (Pirk et
%    al. 2022) with AMIS (Aalstad and Alonso-Gonzalez et al., 2025) as
%    described in Aalstad et al. (2025). 

Np=size(pros,1);
No=size(pred,1);
Ne=size(pros,2);
Nq=size(proc,3);

% Negative log of the target unnormalized posterior
phi=zeros(Ne,Nq);

% Normalizing constant of DM proposals
lc=zeros(Nq,1);
% Log-sum-exp of the "DM" of proposals including normalizing constants
lsepsi=zeros(Ne,Nq); 

isdiag=true;

% Note here ell indexes both assimilation cycles within an IES and
% generations of the IES. So e.g. with Na=4 and Ng=2 then ell=6
% corresponds to the prior of the second generation of IES sine each
% generation will have (Na+1)=5 iterations. 
for ell=1:Nq
    prosell=pros(:,:,ell); % Samples from proposal ell "sampell" (not a typo)
    A0ell=prosell-prim;
    phi0ell=-0.5.*smahal(A0ell,priv,isdiag)...
        -(0.5*Np)*log(2*pi)...
        -0.5.*sum(log(priv));
    predell=pred(:,:,ell);
    resell=obs-predell;
    phidell=-0.5.*smahal(resell,Rvec,isdiag)...
        -(0.5*No)*log(2*pi)...
        -0.5.*sum(log(Rvec));
    phiell=phi0ell+phidell;
    phi(:,ell)=phiell;


    % The DM of proposals denominator term
    psi=zeros(Ne,Nq);
    for j=1:Nq
        mj=prom(:,j);
        Cj=proc(:,:,j);
        Cj=Cj+1e-6.*eye(size(Cj));
        if ell==1 % Compute the proposal normalizing constant once.
            lcj=-(0.5*Np)*log(2*pi)-0.5*logdetc(Cj);
            lc(j)=lcj;
        else % No need to recompute the proposal's normalizing constant.
            lcj=lc(j);
        end
        Aj=prosell-mj;
        psij=-0.5.*smahal(Aj,Cj,false)...
            +lcj-log(Nq);
        psi(:,j)=psij;
    end
    % LSE over all proposals in the mixture evaluated for prosell
    lsepsiell=logsumexp(psi,2);  
    lsepsi(:,ell)=lsepsiell;
end

% Log of unnormalized weights
logwu=phi-lsepsi; % Keep this as Ne x Nq for simplicity?

% Log evidence calculation
Nh=Ne*Nq;
logZ=logsumexp(logwu(:),1)-log(Nh);

% Log of normalized weights
logw=logwu-logZ-log(Nh); 
w=exp(logw);



    function smd=smahal(Ain,Cin,isdiag)
        % Calculates the squared Mahalanobis distance of all the M column
        % vectors (N x 1) contained in Ain (N x M) from the distribution
        % with mean 0 and (N x N) covariance matrix given by Cin.
        if isdiag
            smd=sum((Ain.^2)./Cin,1)';
        else % Recycling chol might be faster but only worth it for large N
            B=Cin\Ain;
            smd=sum(Ain.*B,1)';
        end
    end

    function ldc=logdetc(Cin)
        % Calculates the log of the determinant of a non-diagonal covariance
        % matrix using the Cholesky decomposition.
        % Reuse the Cholesky decomposition for faster inversion in smahal.
        % See e.g. Murphy (2021) PML1 Page 231
        L=chol(Cin,'lower');
        ldc=2*sum(log(diag(L)));
    end

    function lse=logsumexp(a,dim)
        % A simple implementation of the "log-sum-exp" function that is
        % valid for matrices.
        amax=max(a,[],dim);
        ad=a-amax;
        lse=amax+log(sum(exp(ad),dim));
    end

end
