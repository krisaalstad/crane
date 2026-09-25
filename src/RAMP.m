function chain = RAMP(p,f,o,psi0,C0,chaincounter)
% RAMP: Particle Markov Chain Monte Carlo (PMCMC) using the Robust
% Adaptive Metropolis (RAM) algorithm. This implementation is
% adapted to hierarchical crysopheric reanalysis

% Turn off recycle for RAMP, just want accuracy not efficiency
p.recycle=0; 
Nc=1e4; 
Nr=1e2; % Number of posterior samples to subsample from the full chain
burnf=0.5; % Initial fraction to discard as burn-in
Nhy=2*p.Np; % Number of hyperparameter to infer with RAMP
psi=zeros(Nhy,Nc); % psi vector for updating with RAMP
gen.r=[]; % Initialize with empty result struct
gen.p=p;
Nwy=numel(p.wys);

% Doesn't leave this function (default=0).
p.samplingchains=0; % Flag for PMCMC debugging.
p.dostoch=0; % Use DEnKF to memory leakage associated with randn across years

% Could use gradients here too, e.g. Langevin Monte Carlo

% Log prior convenience function
logpri=@(mu,tau,m,s,eta,chi) -(0.5.*sum(((mu-m)./s).^2,1))...
    -(0.5.*sum(((tau-eta)./chi).^2,1));
thesep=1:p.Np;
% m: Hyperprior mean on prior mean mu
hm=p.m(thesep);
% s: Hyperprior std on prior mean mu
hs=p.s(thesep);
% eta: Hyperprior mean on prior scale tau
he=p.e(thesep);
% chi: Hyperprior mean on prior scale tau
hx=p.x(thesep);

% RAM hyperparameters based on Vihola's paper.
crash=0;
crashpsi=[];
mhopt=0.234;
gam=2.0/3.0;%2.0/3.0;
Eye=eye(Nhy);
if isempty(C0)
    sigp=0.3;
    % The initial covariance could also be specified externally
    Cp=(sigp.^2).*Eye; % Initial proposal covariance
else
    Cp=C0;
end
Sc=chol(Cp,'lower');
accepted=0;
curp=0;
arate=accepted/Nc; % Accepting initialization by default

monkeytacos=123; % xkcd reference
rng(chaincounter*monkeytacos); % For reproducibility each chain
% has its own fixed seed.

for c=1:Nc
    percc=round(1e2*c/Nc);
    

    % Check for memory leaks
    if mod(c,100)==0
        [~,vm] = system(sprintf('grep VmRSS /proc/%d/status', feature('getpid')));
        w = whos;
        fprintf('c=%5d  whos=%7.1f MB  %s', c, sum([w.bytes])/2^20, vm);
    end

    if percc>=curp
        fprintf('\n PMCMC complete %d, iteration %d, arate %4.2f, chain %d \n',...
            curp,c,arate,chaincounter)
        curp=curp+1;
    end
    if c==1
        pinit=1; % Perturbed initialization, true or false
        if isempty(psi0)
            % Intialize with the prior if psi0 empty (alt: initialize with e.g. EM)
            psi0=zeros(Nhy,1);
            psi0(1)=p.ma+pinit.*p.sa.*randn(1,1);
            psi0(2)=p.mb+pinit.*p.sb.*randn(1,1);
            psi0(3)=p.mc+pinit.*p.sc.*randn(1,1);
            psi0(p.Np+1)=p.ea+pinit.*p.xa.*randn(1,1);
            psi0(p.Np+2)=p.eb+pinit.*p.xb.*randn(1,1);
            psi0(p.Np+3)=p.ec+pinit.*p.xc.*randn(1,1);
            if ~p.glacier
                psi0(p.Np)=p.mv+pinit.*p.sv.*randn(1,1);
                psi0(Nhy)=p.ev+pinit.*p.xv.*randn(1,1);
            end
        else
            r = randn(Nhy,1);
            psi0=psi0+pinit.*(Sc*r); % Perturb initialization
        end
        h.mua=psi0(1);
        h.mub=psi0(2);
        h.muc=psi0(3);
        h.taua=psi0(p.Np+1);
        h.taub=psi0(p.Np+2);
        h.tauc=psi0(p.Np+3);
        if ~p.glacier
            h.muv=psi0(p.Np);
            h.tauv=psi0(Nhy);
        end

        psi(:,1)=psi0;
        psic=psi0;
    else
        % Propose a move after initializing (c>1)
        r = randn(Nhy,1);
        prop = Sc*r;
        psip = psic+prop;
        logprip=logpri(psip(1:p.Np),psip((p.Np+1):(Nhy)),hm,hs,he,hx);
        h.mua=psip(1);
        h.mub=psip(2);
        h.muc=psip(3);
        h.taua=psip(p.Np+1);
        h.taub=psip(p.Np+2);
        h.tauc=psip(p.Np+3);
        if ~p.glacier
            h.muv=psip(p.Np);
            h.tauv=psip(Nhy);
        end
    end

    
    % Only need to run this for new steps, not on rejection
    try
        cs=rng; % Exctract seed state used for proposal evidence
        gen.r=CRA(p,f,o,h,gen.r,1);
        % Extract log joint evidence (marignal likelihood)
        logevip=gen.r.logZ;
        if c>1
            hc.mua=psic(1);
            hc.mub=psic(2);
            hc.muc=psic(3);
            hc.taua=psic(p.Np+1);
            hc.taub=psic(p.Np+2);
            hc.tauc=psic(p.Np+3);
            if ~p.glacier
                hc.muv=psic(p.Np);
                hc.tauv=psic(Nhy);
            end
            rng(cs); % Reset seed state to same as for the proposal
            gen.r=CRA(p,f,o,hc,gen.r,1); % Rerun current
            logevic=gen.r.logZ;
        else
            logevic=logevip; % Ensures initialization is accepted
        end
    catch ME
        % CRA can crash for certain proposals with unrealistic
        % hyperparameter values, this is a quick fix that can be avoided
        % by initializing at a MAP estimate.
        crash=crash+1;
        fprintf('\n Chain %d crash %d at c=%d: %s (%s line %d)\n', ...
        chaincounter, crash, c, ME.message, ...
        ME.stack(1).name, ME.stack(1).line);
        fprintf('\n Chain %d crashed %d out of %d steps',chaincounter,crash,c);
        crashpsi=[crashpsi psip];
        logevip=-Inf; 
    end
    % Initialize current log of current evidence (marginal likelihood)
    % and the current negative log posterior
    if c==1
        psip=psic; % Ensures acceptance of initialization for convenience
        logprip=logpri(psi0(1:p.Np,c),psi0((p.Np+1):(Nhy),c),hm,hs,he,hx);
        Uc=-(logevip+logprip);
    else
        Uc=-(logevic+logpric);
    end

    Up=-(logevip+logprip); % Proposed negative log posterior
    A=min(1,exp(-Up+Uc));
    u=rand(1);
    accept=A>=u;
    if c==1
        % This is true by construction, but here for clarity
        accept=1; % Always accept initial point
    end

    if accept
        psic=psip;
        logpric=logprip; % This is deterministic
        accepted=accepted+1;
    end
    arate=accepted/c;

    psi(:,c)=psic;


    adapt=1;
    if c>1&&adapt % Adapt if not initial point
        % Always adapt regardless of acceptance or rejection
        stepc=c-1; % Only adapt from c>1
        eta=min(1,Nhy*stepc.^(-gam));
        rinner=(r')*r;
        if A>=mhopt
            % Cholesky update
            Rc=Sc';
            Rc=cholupdate(Rc,sqrt((eta*(A-mhopt))/rinner)*prop);
            Sc=Rc';
        else
            % Cholesky downdate
            Rc=Sc';
            try
                Rc=cholupdate(Rc,sqrt((eta*abs(A-mhopt))/rinner)*prop,'-');
            catch
                Rc=Sc';
                fprintf('\n Downdate failed \n');
            end
            Sc=Rc';
        end
    end

end


% Discard burn in, randomly subsample, and rerun
% CRA for each of these. We do not cross polinate
% with MAGPIES, each psi sample should be run independently
% through CRA.

burn=round(burnf*Nc);
psiburn=psi(:,1:burn);
psi=psi(:,(burn+1):Nc);
psifull=psi;

% Redefine chain length after discarding burn in
Nc=size(psi,2);
Nr=min(Nr,Nc);
if Nr<Nc
    resample=randsample(Nc,Nr,false);
    psi=psi(:,resample); % These are the final samples for postpred
end


chain.psi=psi;
chain.psifull=psifull;
chain.psiburn=psiburn;
% The burn in can be "rescued" later if needed
chain.arate=arate;


p.dostoch=1; % Reset dostoch


end
