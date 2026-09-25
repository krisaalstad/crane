%% Create real experiment for a particular site.
clearvars;
% Add src to the path if needed.
if ~any(which('MAGPIES'))
    paths;
end
setpars;
close all;


p.sigD=1e1; % SWE obs error std (mm)
p.sigF=0.15; 
p.dostoch=1;


p.vispri=1;
p.lowy=1; % Explicitly leave out data to create a validation period
p.doMCMC=0; % Benchmarking. Warning! This will slow down the script (a lot).
p.assimD=1; % Do SWE assimilation
p.assimF=0; % Do FSCA assimilation

% Check input data exists, otherwise fetch it from zenodo
if ~any(exist('input','dir'))
    get_input;
end

if p.doMCMC
    Nchains=10; %=number of cores, make sure you have these.
    % Tested  on Matlab R2021a
    p.dopar=1; 
    fprintf('\n Running new PMCMC with %d chains parallel=%d',Nchains,p.dopar);
end

p.recycle=1;
p.Ng=2;
p.Ne=100; 
p.Na=4;



snowsites={'HS'; 'AK'; 'FF'; 'WF'};

for ss=1:4

rng(54321+ss); % For reproducibility
snowsite=snowsites{ss};

fprintf('\n Running %s assimF%d \n',snowsite,p.assimF);

obsfile=sprintf('input/snow/%s/%s_obs.mat',snowsite,snowsite);
load(obsfile);
ffile=sprintf('input/snow/%s/%s_forcing.mat',snowsite,snowsite);
load(ffile);



% Specify the period, should correspond to water years
otv=datevec(obs.t);
wystart=otv(1,1)+1; 
wyend=max(otv(:,1));
wys=wystart:wyend; wys=wys';
tstart=datenum(sprintf('01-Sep-%d',wystart-1));
tend=datenum(sprintf('01-Sep-%d',wyend))-1;
t=tstart:tend; t=t'; Nt=numel(t);


if p.assimD&&p.lowy
    p.lowys=(2000:wyend)'; % Left out water years (validation period).
else % For FSCA the calibraiton is always just 2000-2023
    p.lowy=0; % Can also be off for SWE, then all data are assimilated.
    p.lowys=[];
end


theseo=obs.t>=tstart&obs.t<=tend;
obs.D=obs.D(theseo); obs.t=obs.t(theseo);
% Assimilate ~fortnightly (10 days at higher latitude) SWE obs
if ss<4 % WF already fortnightly 
    theset=zeros(numel(obs.D),'logical');
    told=-Inf;
    for j=1:numel(obs.t)
        tj=obs.t(j);
        dt=tj-told;
        if dt>=10
            theset(j)=1;
            told=tj;
        end
    end
    obs.D=obs.D(theset);
    obs.t=obs.t(theset);
end
theseos=obs.ts>=tstart&obs.ts<=tend;
obs.Fs=obs.Fs(theseos);
obs.ts=obs.ts(theseos);
obs.vzs=obs.vzs(theseos);
otvs=datevec(obs.ts);
motvs=otvs(:,2);
here=motvs>=0&motvs<13; % Optional month filter, not used.
obs.Fs=obs.Fs(here);
obs.ts=obs.ts(here);
obs.vzs=obs.vzs(here);
% Remove high view zenith (>60^deg)
highz=obs.vzs>60; % 
obs.Fs=obs.Fs(~highz);
obs.ts=obs.ts(~highz);


thesef=f.t>=tstart&f.t<=tend;
f.P=f.P(thesef);
f.T=f.T(thesef);
f.t=f.t(thesef);

p.wys=wys;
Nwy=numel(wys);


skip=0; % Skip to PMCMC, just for testing
if ~skip


%% Particle Expectation Maximization    
tic;
gen = particleEM(p,f,obs);
ttoc=toc;
fprintf('\n Ran %d generation particleEM in %4.2f seconds \n',p.Ng,ttoc);


%% Laplace approximation
Q=gen.Hhat; % Precision matrix = Hessian
% Check that the Laplace approximation behaves well numerically with no
% non positive definite Hessian issues. 
try
    L=chol(Q,'lower');
    CLap=inv(Q);
    pdhess=1;
catch
    error('NPD Hessian');
    hpris=[p.s(1:4); p.x(1:4)]; 
    Q=inv(diag(hpris.^2));
    L=chol(Q);
    CLap=diag(hpris.^2);
    pdhess=0;
end


%% Perform AdaPBS hyperparameter inference using the WAMIS function

% Note that the EM is not just used for the Laplace approximation, that is
% just a by-product. The main purpose is to guide MAGPIES.
Ns=1e2; % Sample 100 hyperparameter particles
% Start by sampling from a mixture of the Laplace approximation and the
% hyperprior.
Nhyp=p.Np*2;
psim=gen.psi(:,end,end);
psip=psim+(L')\randn(size(L,1),Ns);
% Here we also effectively use a deterministic mixture, sampling equally many form
% each of the "subproposals" i.e. 100 from prior and 100 from Laplace Approx
psiLaplace=psim+(L')\randn(size(L,1),Ns);
hprim=[p.m(1:p.Np); p.e(1:p.Np)];
hpris=[p.s(1:p.Np); p.x(1:p.Np)]; 
hpriv=hpris.^2;
psihpri=hprim+hpris.*randn(Nhyp,Ns);
diversity=50;
Nr=50;
psi=zeros(Nhyp,Ns,Nr);
hprom=zeros(Nhyp,Nr);
hprom(:,1)=hprim;
hprom(:,2)=psim;
hproc=zeros(Nhyp,Nhyp,Nr);
hproc(:,:,1)=diag(hpriv);
hproc(:,:,2)=CLap; 
logZ=zeros(Ns,Nr);
adapt=0.1;
tic;
r=0;
doadapt=1;
while doadapt
    r=r+1;
    if r==1
        psi(:,:,r)=psihpri;
    elseif r==2
            psi(:,:,r)=psiLaplace;
    end
    % Note, we only need to run reMAGPIES for the newly proposed psi(:,:,r)
    [~,lZn] = reMAGPIES(p,psi(:,:,r),gen);
    logZ(:,r)=sum(lZn,1);

    % Only need WAMIS for r>1 (prior), use Laplace proposal for r=2
    if r>1

        % WAMIS update
        [omega,~] = WAMIS( logZ(:,1:r), hprim, hpriv,...
            psi(:,:,1:r), hprom(:,1:r), hproc(:,:,1:r) );
        psihist=psi(:,:,1:r);
        omega=omega(:);
        Neff=1./sum(omega.^2);

        % Early stop if Neff>=diversity
        if Neff>diversity
            converged=1;
            disp('converged early!')
            Nr=r;
        end

        doadapt=r<Nr; 
        Nh=r*Ns; % Total number of historical samples from current DM proposal

        if doadapt
            clip=round(adapt*Ns);
            omegas=sort(omega,'descend');
            omegac=omegas(clip);
            if omegac>0
                omega(omega>omegac)=omegac;
                omega=omega/sum(omega);
            end
        end

        psihist=reshape(psihist,[Nhyp,Nh]);
        resample=randsample(Nh,Ns,true,omega);
        psire=psihist(:,resample);

        if doadapt
            pm=mean(psire,2);
            hprom(:,r+1)=pm;
            Z=randn(Nhyp,Ns);
            if omegac>0 % Adapt covariance if not degenerate 
                A=psire-pm;
                pc=(1/Ns).*(A*A');
            else % Otherwise use prior covariance
                fprintf('\n Degenerate \n')
                pc=diag(hpriv).^(1/r);
            end
            hproc(:,:,r+1)=pc;
            try
                L=chol(pc,'lower');
            catch
                pc=pc+1e-3.*eye(size(pc));
                L=chol(pc,'lower');
            end
            psip=pm+L*Z;
            psi(:,:,r+1)=psip;
        end
        
    end
end
ttoc=toc;
fprintf('\n %d WAMIS iterations took %4.2f seconds Neff=%d \n',...
    r,ttoc,round(Neff));

%% Posterior predictive simulations.
resample=randsample(Nh,p.Ne,true,omega);

% Speed up trick: Just reweight the unique resamples by their frequency
[c,~,ic]=unique(resample);
counts=accumarray(ic,1);
uweight=counts/sum(counts);
uresample=c;
psipostu=psihist(:,uresample);
Nu=numel(uweight);
Dpo=zeros(Nt,p.Nr,Nu);
mweight=uweight./p.Nr; % Merged weights (inner are equal)

pp=p;
pp.case='PP';
pp.recycle=0;
pp.resample=1;
pp.Na=4;
pp.Ng=1;
pp.Ne=1e2;

tic;
popred=postpred(pp,f,obs,psipostu,mweight);
ttoc=toc;
fprintf('\n Posterior predictions took %4.2f seconds \n',ttoc);
Dpom=popred.Dm;
Dpos=popred.Ds;

if p.assimF
    Fpom=popred.Fm;
    Fpos=popred.Fs;
else
    Fpom=[];
    Fpos=[];
end

hpri=1;
if hpri==1
    Ndraw=Nu;
else
    Ndraw=1;
end
psipri=zeros(2*p.Np,Ndraw);
perturb=(Ndraw>1);
psipri(1,:)=p.ma+perturb.*p.sa.*randn(1,Ndraw);
psipri(2,:)=p.mb+perturb.*p.sb.*randn(1,Ndraw);
psipri(3,:)=p.mc+perturb.*p.sc.*randn(1,Ndraw);
psipri(4,:)=p.mv+perturb.*p.sv.*randn(1,Ndraw);
psipri(5,:)=p.ea+perturb.*p.xa.*randn(1,Ndraw);
psipri(6,:)=p.eb+perturb.*p.xb.*randn(1,Ndraw);
psipri(7,:)=p.ec+perturb.*p.xc.*randn(1,Ndraw);
psipri(8,:)=p.ev+perturb.*p.xv.*randn(1,Ndraw);

pp.Na=0; % Ensures only prior simulations
tic;
prpred=postpred(pp,f,obs,psipri,[]);
ttoc=toc;
fprintf('\n Prior predictions took %4.2f seconds \n',ttoc);
Dprm=prpred.Dm;
Dprs=prpred.Ds;
if p.assimF
    Fprm=prpred.Fm;
    Fprs=prpred.Fs;
else
    Fprm=[];
    Fprs=[];
end

% Create a "pppred" struct
pppred=popred;
pppred.Dprim=Dprm;
pppred.Dpris=Dprs;
if p.assimF
    pppred.Fprim=Fprm;
    pppred.Fpris=Fprs;
end

pppred.psipri=psipri;
pppred.psi=psipostu;
pppred.weight=uweight;

set(0, 'DefaultFigureRenderer', 'painters');
flag=snowpercplot(t,pp,obs,Dprm,Dprs,Dpom,Dpos,...
    Fprm,Fprs,Fpom,Fpos,snowsite);

%% Try CP

pcp=p;
pcp.case='CP';
pcp.resample=1;
pcp.Ne=1e3;
pcp.Ng=1;
pcp.Nr=1e2;
pcp.Na=4;
psicp=zeros(2*p.Np,1);
psicp(1)=p.ma;
psicp(2)=p.mb;
psicp(3)=p.mc;
psicp(4)=p.mv;
psicp(5)=p.ea;
psicp(6)=p.eb;
psicp(7)=p.ec;
psicp(8)=p.ev;
tic;
cppred=postpred(pcp,f,obs,psicp,[]);
ttoc=toc;
fprintf('\n CP predictions took %4.2f seconds \n',ttoc);
Dpom=cppred.Dm;
Dpos=cppred.Ds;
if p.assimF
    Fpom=cppred.Fm;
    Fpos=cppred.Fs;
else
    Fpom=[];
    Fpos=[];
end

Dprm=cppred.Dprim;
Dprs=cppred.Dpris;
if p.assimF
    Fprm=cppred.Fprim;
    Fprs=cppred.Fpris;
else
    Fprm=[];
    Fprs=[];
end

cppred.psi=psicp;



flag=snowpercplot(t,pcp,obs,Dprm,Dprs,Dpom,Dpos,...
    Fprm,Fprs,Fpom,Fpos,snowsite);

%% Try NP

p.vispri=1;
pnp=p;
pnp.case='NP';
pnp.Ng=1;
pnp.recycle=0;
pnp.resample=1;
pnp.Ne=1e3;
pnp.Nr=1e2;
pnp.Na=4;
psinp=zeros(2*p.Np,1);
psinp(1)=p.ma;
psinp(2)=p.mb;
psinp(3)=p.mc;
psinp(4)=p.mv;
psinp(5)=p.ea;
psinp(6)=p.eb;
psinp(7)=p.ec;
psinp(8)=p.ev;
tic;
nppred=postpred(pnp,f,obs,psinp,[]);
ttoc=toc;
fprintf('\n NP predictions took %4.2f seconds \n',ttoc);
Dpom=nppred.Dm;
Dpos=nppred.Ds;
if p.assimF
    Fpom=nppred.Fm;
    Fpos=nppred.Fs;
else
    Fpom=[];
    Fpos=[];
end

Dprm=nppred.Dprim;
Dprs=nppred.Dpris;
if p.assimF
    Fprm=nppred.Fprim;
    Fprs=nppred.Fpris;
else
    Fprm=[];
    Fprs=[];
end

nppred.psi=psinp;

flag=snowpercplot(t,pnp,obs,Dprm,Dprs,Dpom,Dpos,...
    Fprm,Fprs,Fpom,Fpos,snowsite);

%% Try HMAP (NP with optimized hyperpars)
pmap=p;
pmap.vispri=1;
pmap.case='MAP';
pmap.Ng=1;
pmap.Na=4;
pmap.recycle=0;
pmap.resample=1;
pmap.Ne=1e2;
pmap.Nr=1e2;
psimap=squeeze(gen.psi(:,end,end));
tic;
mappred=postpred(pmap,f,obs,psimap,[]);
ttoc=toc;
fprintf('\n MAP predictions took %4.2f seconds \n',ttoc);

Dpom=mappred.Dm;
Dpos=mappred.Ds;
if p.assimF
    Fpom=mappred.Fm;
    Fpos=mappred.Fs;
else
    Fpom=[];
    Fpos=[];
end


% Assign real prior (not tuned) to map from prpred
Dprm=prpred.Dm;
Dprs=prpred.Ds;
mappred.Dprim=Dprm;
mappred.Dpris=Dprs;
if p.assimF
    Fprm=prpred.Fm;
    Fprs=prpred.Fs;
    mappred.Fprim=Fprm;
    mappred.Fpris=Fprs;
else
    Fprm=[];
    Fprs=[];
end
mappred.psi=psimap;
mappred.CLap=CLap;

flag=snowpercplot(t,pmap,obs,Dprm,Dprs,Dpom,Dpos,...
    Fprm,Fprs,Fpom,Fpos,snowsite);

%% Save the results

clear r;
r.map=mappred;
r.cp=cppred;
r.np=nppred;
r.pp=pppred;
r.p=p;
r.obs=obs;

tofile=sprintf('results/Snow_%s_D%d_F%d_LO%d.mat',...
    snowsite,p.assimD,p.assimF,p.lowy);

save(tofile,'r');
end

%% PMCMC

if p.doMCMC
    pmcmc=p;
    pmcmc.case='PP';
    pmcmc.recycle=0;
    pmcmc.resample=1;
    pmcmc.Na=4;
    pmcmc.Ng=1;
    pmcmc.Ne=1e2;
    dospeedup=1;
    if skip
        dospeedup=0;
    end
if dospeedup
    psi0=psimap;
    if pdhess==1
        C0=CLap;
    else
        C0=[];
    end
else
    psi0=[];
    C0=[];
end

% Clean up workspace
if ~skip
    clear Dpo gen wn
end

pname='Processes'; 
par=gcp('nocreate');

% Requires older (R2021a) Matlab...
tic;
if Nchains>1
    tarfile=sprintf('results/chains_%s_L%d_D%d_F%d.mat',...
        snowsite,p.lowy,p.assimD,p.assimF);
    if ~any(exist(tarfile,'file'))
        if isempty(par)
            if p.dopar
                delete(par);
                par=parpool(pname,Nchains);
            end
        elseif ~par.Connected||par.NumWorkers~=Nchains&&p.dopar
            delete(par);
            par=parpool(pname,Nchains);
        end
        if p.dopar
            parfor j=1:Nchains
                chains(j).chain=RAMP(pmcmc,f,obs,psi0,C0,j);
            end
        else
            for j=1:Nchains
                chains(j).chain=RAMP(pmcmc,f,obs,psi0,C0,j);
            end
        end
        save(tarfile,'chains');
        delete(par);
    else
        fprintf('\n Loading chains existing file \n')
        load(tarfile);
    end
    chain.psi=[];
    for j=1:Nchains
        chain.psi=[chain.psi chains(j).chain.psifull];
    end
    Ncc=size(chain.psi,2); % Size of concatenated chain
    if p.Nss<Ncc
        fprintf('\n Resampling chains \n');
        resample=randsample(Ncc,p.Nss,false);
        chain.psi=chain.psi(:,resample); % These are the final samples for postpred
    end
else
    chain = RAMP(pmcmc,f,obs,psi0,C0,j); 
end
ttoc;
fprintf('\n PMCMC took %4.2f seconds \n',ttoc);

psipmcmc=chain.psi;


% Here outer weights are equal and inner weights are equal
mweight=1/(size(psipmcmc,2)*p.Nr);

tic;
pmcmc.Na=4;
pmcmc.dostoch=0;
pmcmcpred=postpred(pmcmc,f,obs,psipmcmc,mweight);
ttoc=toc;
fprintf('\n Posterior PMCMC predictions took %4.2f seconds \n',ttoc);
Dpom=pmcmcpred.Dm;
Dpos=pmcmcpred.Ds;

if p.assimF
    Fpom=pmcmcpred.Fm;
    Fpos=pmcmcpred.Fs;
else
    Fpom=[];
    Fpos=[];
end
pmcmcpred.psi=psipmcmc;
pmcmcpred.pri='Same prior as PP'; 
pmcmcpred.chain=chain;

Dprm=prpred.Dm;
Dprs=prpred.Ds;
pmcmcpred.Dprm=Dprm;
pmcmcpred.Dprs=Dprs;
if p.assimF
    Fprm=prpred.Fm;
    Fprs=prpred.Fs;
    pmcmcpred.Fprm=Fprm;
    pmcmcpred.Fprs=Fprs;
else
    Fprm=[];
    Fprs=[];
end

set(0, 'DefaultFigureRenderer', 'painters');
pmcmc.case='PMCMC'; % For plotting
flag=snowpercplot(t,pmcmc,obs,Dprm,Dprs,Dpom,Dpos,...
    Fprm,Fprs,Fpom,Fpos,snowsite);

r.mcmc=pmcmcpred;
save(tofile,'r');

end





end
