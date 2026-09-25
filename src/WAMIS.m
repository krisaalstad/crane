function [w,logZ] = WAMIS( plml, hprim, hpriv, hpros, hprom, hproc )
%% WAMIS w = WAMIS( pll, hprim, hpriv, hpros, hprom, hproc )
%
% NB: This algorithm is the same as AMIS but uses a particle approximation
%     of the marginal likelihood that is provided as an input. As such the
%     WAMIS scheme is ideally suited as a kind of pseudo-marginal method 
%     for inferring hyperparameters in hierarchical infernece via nested
%     particle methods.
%
% Dimensions: 
%             Np = Number of HYPERparameters to update.
%             Ne = Number of ensemble members (particles) per iteration
%             Nq = Number of proposals om the det. mixture = (Na+1) x Ng
%             Nh = Number of historical particles = Ne x Nq
% 
% Inputs:
%    plml: Particle log marginal likelihood approximation (Ne x Nq)
%    hprim: The mean vector of the hyperprior (Np x 1 array)
%    hpriv: The variance vector of the hyperprior (Np x 1 array)
%    hpros: Samples from the proposal (Np x Ne x Nq arrray)
%    hprom: Mean of the proposal (Np x Nq array)
%    hproc: Covariance of the proposal (Np x Np x Nq array)
%        
%    Outputs:
%        w: Posterior weights (Nh x 1 array)
%        logZ : Log evidence (marginal likelihood)
%    Dimensions

Np=size(hpros,1);
Ne=size(hpros,2);
Nq=size(hproc,3);

% Negative log of the target unnormalized hyperposterior
phi=zeros(Ne,Nq);

% Normalizing constant of DM proposals
lc=zeros(Nq,1);
% Log-sum-exp of the "DM" of proposals including normalizing constants
lsepsi=zeros(Ne,Nq); 

for ell=1:Nq
    prosell=hpros(:,:,ell); % Samples from proposal ell 
    A0ell=prosell-hprim;
    phi0ell=-0.5.*smahal(A0ell,hpriv,true)...
        -(0.5*Np)*log(2*pi)...
        -0.5.*sum(log(hpriv));
    phidell=plml(:,ell);
    phiell=phi0ell+phidell;
    phi(:,ell)=phiell;


    % The DM of proposals denominator term
    psi=zeros(Ne,Nq);
    for j=1:Nq
        mj=hprom(:,j);
        Cj=hproc(:,:,j);
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
            B=Ain./Cin;
        else % Recycling chol might be faster but only worth it for large N
            B=Cin\Ain;
        end
        smd=sum(Ain.*B,1)';
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
