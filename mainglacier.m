%% Create real experiment for a particular site.
clearvars;
glacier=1;
setpars;
close all;

% Custom parameters for this experiment
p.sigD=2e2; % Obs. error std for winter MB, annual MB is double this.
p.doMCMC=0; % Benchmarking switch. Warning! This will slow down the script (a lot).
p.visD=0;
p.lowy=1; % Explicitly leave out data to create a validation period


p.recycle=1;
p.dostoch=1;
p.Ng=2;
p.Na=4;
p.Ne=1e2;

% Check input data exists, otherwise fetch it from zenodo
if ~any(exist('input','dir'))
    get_input;
end

if p.doMCMC
    Nchains=10; % =ncores, make sure you have these.
    % Parallel chains only tested on Matlab R2021a.
    p.dopar=1;
    fprintf('\n Running PMCMC with %d chains parallel=%d',Nchains,p.dopar);
end

if p.lowy
    p.lowys=(2000:2024)'; % Left out water years (validation period)
else
    p.lowys=[];
end

glaciers={'AB'; 'SG'; 'SB'; 'CD'};

for gg=1:4


rng(12345+gg); % For reproducibility
glacier=glaciers{gg};
fprintf('\n Running %s \n',glacier);

obsfile=sprintf('input/glacier/%s/%s_obs.mat',glacier,glacier);
load(obsfile);
ffile=sprintf('input/glacier/%s/%s_forcing.mat',glacier,glacier);
load(ffile);

obs.D=obs.B.*1e3; % to mm

% Specify the period, should correspond to water years
otv=datevec(obs.t);
wystart=otv(1,1);
wyend=max(otv(:,1));
wys=wystart:wyend; wys=wys';
tstart=datenum(sprintf('01-Sep-%d',wystart-1));
tend=datenum(sprintf('01-Sep-%d',wyend))-1;
t=tstart:tend; t=t'; Nt=numel(t);

% Homogenize the annual MB obs to all be end of August for simplicity.
% Can later just use varying length DAWs (should work with gen structure)
% for glaciers. This is to properly capture the accumulation from the annual SMB
% date (summer balance) to the subsequent winter balance date. The
% homogenization approach here is a quick fix that ensures at least the
% modeled cumulative mass balance will be correct.
otv(2:2:end,2)=8;
otv(2:2:end,3)=31;
obs.t=datenum(otv);
theseo=obs.t>=tstart&obs.t<=tend;
obs.D=obs.D(theseo); obs.t=obs.t(theseo);

thesef=f.t>=tstart&f.t<=tend;
f.P=f.P(thesef);
f.T=f.T(thesef);
f.t=f.t(thesef);

p.wys=wys;
p.assimD=1;
p.assimF=0;
Nwy=numel(wys);


%% Run particle Expectation Maximization
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
    hpris=[p.s(1:3); p.x(1:3)]; 
    Q=inv(diag(hpris.^2));
    CLap=diag(hpris.^2);
    L=chol(Q,'lower');
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
    [wn,lZn] = reMAGPIES(p,psi(:,:,r),gen);
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
        if Neff>=diversity
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
Cpom=popred.Cm;
Cpos=popred.Cs;
pppred=popred;

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
psipri(4,:)=p.ea+perturb.*p.xa.*randn(1,Ndraw);
psipri(5,:)=p.eb+perturb.*p.xb.*randn(1,Ndraw);
psipri(6,:)=p.ec+perturb.*p.xc.*randn(1,Ndraw);

pp.Na=0; % Ensures only prior simulations
tic;
prpred=postpred(pp,f,obs,psipri,[]);
ttoc=toc;
fprintf('\n Prior predictions took %4.2f seconds \n',ttoc);
Dprm=prpred.Dm;
Dprs=prpred.Ds;
Cprm=prpred.Cm;
Cprs=prpred.Cs;

pppred.Dprim=Dprm;
pppred.Dpris=Dprs;
pppred.Cprim=Cprm;
pppred.Cpris=Cprs;

pppred.psipri=psipri;
pppred.psi=psipostu;
pppred.weight=uweight;



set(0, 'DefaultFigureRenderer', 'painters');
flag=glacierpercplot(t,pp,obs,Dprm,Dprs,Dpom,Dpos,...
    Cprm,Cprs,Cpom,Cpos,glacier);

%% Try CP
pcp=p;
pcp.case='CP';
pcp.resample=1;
pcp.Ne=1e3;
pcp.Na=4;
pcp.Ng=1;
pcp.Nr=1e2;
psicp=zeros(2*p.Np,1);
psicp(1)=p.ma;
psicp(2)=p.mb;
psicp(3)=p.mc;
psicp(4)=p.ea;
psicp(5)=p.eb;
psicp(6)=p.ec;
tic;
cppred=postpred(pcp,f,obs,psicp,[]);
ttoc=toc;
fprintf('\n CP predictions took %4.2f seconds \n',ttoc);
Dpom=cppred.Dm;
Dpos=cppred.Ds;
Cpom=cppred.Cm;
Cpos=cppred.Cs;
cppred.psi=psicp;

Dprm=cppred.Dprim;
Dprs=cppred.Dpris;
Cprm=cppred.Cprim;
Cprs=cppred.Cpris;

flag=glacierpercplot(t,pcp,obs,Dprm,Dprs,Dpom,Dpos,...
    Cprm,Cprs,Cpom,Cpos,glacier);

%% Try NP

pnp=p;
pnp.case='NP';
pnp.Ng=1;
pnp.Na=4;
pnp.recycle=0;
pnp.resample=1;
pnp.Ne=1e3;
pnp.Nr=1e2;
psinp=zeros(2*p.Np,1);
psinp(1)=p.ma;
psinp(2)=p.mb;
psinp(3)=p.mc;
psinp(4)=p.ea;
psinp(5)=p.eb;
psinp(6)=p.ec;
tic;
nppred=postpred(pnp,f,obs,psinp,[]);
ttoc=toc;
fprintf('\n NP predictions took %4.2f seconds \n',ttoc);
Dpom=nppred.Dm;
Dpos=nppred.Ds;
Cpom=nppred.Cm;
Cpos=nppred.Cs;

Dprm=nppred.Dprim;
Dprs=nppred.Dpris;
Cprm=nppred.Cprim;
Cprs=nppred.Cpris;

nppred.psi=psinp;


flag=glacierpercplot(t,pnp,obs,Dprm,Dprs,Dpom,Dpos,...
    Cprm,Cprs,Cpom,Cpos,glacier);


%% Try HMAP (NP with optimized hyperpars)
pmap=p;
pmap.case='MAP';
pmap.Ng=1;
pmap.Na=4;
pmap.recycle=0;
pmap.resample=1;
pmap.Ne=1e2;
pmap.Nr=1e2;
psimap=squeeze(gen.psi(:,end,end));
mappred=postpred(pmap,f,obs,psimap,[]);
Dpom=mappred.Dm;
Dpos=mappred.Ds;
Cpom=mappred.Cm;
Cpos=mappred.Cs;
mappred.psi=psimap;
Dprm=pppred.Dprim;
Dprs=pppred.Dpris;
Cprm=pppred.Cprim;
Cprs=pppred.Cpris;
mappred.Dprim=Dprm;
mappred.Dpris=Dprs;
mappred.Cprim=Cprm;
mappred.Cpris=Cprs;

tic;
flag=glacierpercplot(t,pmap,obs,Dprm,Dprs,Dpom,Dpos,...
    Cprm,Cprs,Cpom,Cpos,glacier);
ttoc=toc;
fprintf('\n MAP predictions took %4.2f seconds \n',ttoc);

%% Save the results

clear r;


r.map=mappred;
r.cp=cppred;
r.np=nppred;
r.pp=pppred;
r.p=p;
r.obs=obs;

tofile=sprintf('results/Glacier_%s_LO%d.mat',...
    glacier,p.lowy);

save(tofile,'r');

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
clear Dpo gen wn


pname='Processes'; 
par=gcp('nocreate');

tic;
if Nchains>1
    tarfile=sprintf('results/chains_%s_L%d_I%d.mat',...
        glacier,p.lowy,dospeedup);
    if ~any(exist(tarfile,'file'))
        if isempty(par)
            delete(par);
            par=parpool(pname,Nchains);
        elseif ~par.Connected||par.NumWorkers~=Nchains
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
        load(tarfile);
    end
    chain.psi=[];
    for j=1:Nchains
        chain.psi=[chain.psi chains(j).chain.psifull];
    end
    Ncc=size(chain.psi,2); 
    if p.Nss<Ncc
        resample=randsample(Ncc,p.Nss,false);
        chain.psi=chain.psi(:,resample); 
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
pmcmc.dostoch=0;
pmcmc.Na=4;
pmcmcpred=postpred(pmcmc,f,obs,psipmcmc,mweight);
ttoc=toc;
fprintf('\n Posterior PMCMC predictions took %4.2f seconds \n',ttoc);
Dpom=pmcmcpred.Dm;
Dpos=pmcmcpred.Ds;
Cpom=pmcmcpred.Cm;
Cpos=pmcmcpred.Cs;

pmcmcpred.psi=psipmcmc;
pmcmcpred.pri='Same prior as PP'; 
pmcmcpred.chain=chain;

Dprm=prpred.Dm;
Dprs=prpred.Ds;
Cprm=prpred.Cm;
Cprs=prpred.Cs;
pmcmcpred.Dprm=Dprm;
pmcmcpred.Dprs=Dprs;
pmcmcpred.Cprm=Cprm;
pmcmcpred.Cprs=Cprs;

set(0, 'DefaultFigureRenderer', 'painters');
pmcmc.case='PMCMC'; % For plotting
flag=glacierpercplot(t,pmcmc,obs,Dprm,Dprs,Dpom,Dpos,...
    Cprm,Cprs,Cpom,Cpos,glacier);

r.mcmc=pmcmcpred;
save(tofile,'r');

end







end
