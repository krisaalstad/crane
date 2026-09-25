function pred=postpred(p,f,obs,psi,omega)
% Assumes p.recycle=0 and p.resample=1
% if p.Na=0 this will just return the prior predictions


Nu=size(psi,2); % Number of psi vectors to consider 
wys=p.wys;
Nwy=numel(wys);
nothier=any(strcmp(p.case,'CP'))||any(strcmp(p.case,'NP'))||any(strcmp(p.case,'MAP'));
% Include prior prediction for nonhierarchical methods
% (Note that it is also possible to run NP with the MAP-II psi in this way)

logit=@(xb,bd) log(((xb-bd(1))./(bd(2)-bd(1)))...
    ./(1-((xb-bd(1))./(bd(2)-bd(1)))));
expit=@(yu,bd) bd(1)+(bd(2)-bd(1))./(1+exp(-yu));
bds=[p.ab; [-Inf Inf]; p.cb; p.vb];



if p.Na==0||nothier
    omega=1/(p.Nr*Nu);
end


for j=1:Nu
    g=1;
    r=[];
    hyper.mua=psi(1,j);
    hyper.mub=psi(2,j);
    hyper.muc=psi(3,j);
    hyper.taua=psi(p.Np+1,j);
    hyper.taub=psi(p.Np+2,j);
    hyper.tauc=psi(p.Np+3,j);
    if ~p.glacier
        hyper.muv=psi(p.Np,j);
        hyper.tauv=psi(2*p.Np,j);
    end

    r=CRA(p,f,obs,hyper,r,g);
    if j==1
        Nt=size(r.Dpri,1);
        D=zeros(Nt,p.Nr,Nu);
        if any(strcmp(p.case,'CP'))
            phi=zeros(p.Np,p.Nr,Nu);
            phi0=zeros(p.Np,p.Nr);
        else
            phi=zeros(Nwy,p.Np,p.Nr,Nu);
            phi0=zeros(Nwy,p.Np,p.Nr);
        end
        if ~nothier
            phi0=zeros(Nwy,p.Np,p.Nr);
        end
        if ~p.glacier
            F=zeros(Nt,p.Nr,Nu);
        end
    end

    if p.Na>0
        D(:,:,j)=r.Dpos;
        thetaj=r.theta;
        if ~p.glacier
            F(:,:,j)=r.Fpos;
        end
    else
        D(:,:,j)=r.Dpri;
        thetaj=r.theta0;
        if ~p.glacier
            F(:,:,j)=r.Fpri;
        end
    end
    theta0j=r.theta0;
    for k=1:p.Np
        if k~=2
            if any(strcmp(p.case,'CP'))
                phi(k,:,j)=expit(thetaj(k,:),bds(k,:));
                if ~nothier
                    phi0(k,:,j)=expit(theta0j(k,:),bds(k,:));
                end
            else
                phi(:,k,:,j)=expit(thetaj(:,k,:),bds(k,:));
                if ~nothier
                    phi0(:,k,:,j)=expit(theta0j(:,k,:),bds(k,:));
                end
            end
        else
            if any(strcmp(p.case,'CP'))
                phi(k,:,j)=thetaj(k,:);
                if ~nothier
                    phi0(k,:,j)=theta0j(k,:);
                end
            else
                phi(:,k,:,j)=thetaj(:,k,:);
                if ~nothier
                    phi0(:,k,:,j)=theta0j(:,k,:);
                end
            end
        end
    end
end
% Make a function for this to limit memory use
D=permute(D,[3 2 1]);
Dpm=squeeze(sum(sum(omega.*D,1),2));
Dpmr=repmat(Dpm,[1,p.Nr,Nu]);
Dpmr=permute(Dpmr,[3 2 1]);
Dps=sqrt(squeeze(sum(sum(omega.*(D-Dpmr).^2,1),2)));
pred.Dm=Dpm;
pred.Ds=Dps;
D=permute(D,[3 2 1]); % and back (for C)

if any(strcmp(p.case,'CP'))
    phi=permute(phi,[3 2 1]);
    phim=squeeze(sum(sum(omega.*phi,1),2));
    phimr=reshape(phim,1,1,[]);
    phis=sqrt(squeeze(sum(sum(omega.*(phi-phimr).^2,1),2)));
    pred.phim=phim;
    pred.phis=phis;

else
    phi=permute(phi,[4 3 1 2]);
    phim=squeeze(sum(sum(omega.*phi,1),2));
    phimr=repmat(phim,[1,1,p.Nr,Nu]);
    phimr=permute(phimr,[4 3 1 2]);
    phis=sqrt(squeeze(sum(sum(omega.*(phi-phimr).^2,1),2)));
    pred.phim=phim;
    pred.phis=phis;
end


if nothier
    pred.Dprim=mean(r.Dpri,2);
    pred.Dpris=std(r.Dpri,1,2);
    pred.phiprm=mean(phi0,3);
    pred.phiprs=std(phi0,1,3);
end
if ~p.glacier
    F=permute(F,[3 2 1]);
    Fpm=squeeze(sum(sum(omega.*F,1),2));
    Fpmr=repmat(Fpm,[1,p.Nr,Nu]);
    Fpmr=permute(Fpmr,[3 2 1]);
    Fps=sqrt(squeeze(sum(sum(omega.*(F-Fpmr).^2,1),2)));
    pred.Fm=Fpm;
    pred.Fs=Fps;
    if nothier
        pred.Fprim=mean(r.Fpri,2);
        pred.Fpris=std(r.Fpri,1,2);
    end
end

if p.glacier
    C=zeros(Nt,p.Nr,Nu);
    Clast=0;  
    if nothier
        Cpr=zeros(Nt,p.Nr);
        Cprlast=0;
    end
    for n=1:Nwy
        wyk=wys(n);
        tstart=datenum(sprintf('%s-%d',p.dstart,wyk-1));
        tend=datenum(sprintf('%s-%d',p.dend,wyk));
        here=r.t>=tstart&r.t<=tend;
        Dn=D(here,:,:);
        Cn=Clast+Dn;
        C(here,:,:)=Cn;
        Clast=Cn(end,:,:);
        if nothier
            Dprn=r.Dpri(here,:);
            Cprn=Cprlast+Dprn;
            Cpr(here,:)=Cprn;
            Cprlast=Cprn(end,:);
        end
    end
    C=permute(C,[3 2 1]);
    Cpm=squeeze(sum(sum(C,1),2)./(p.Nr*Nu));
    Cpmr=repmat(Cpm,[1,p.Nr,Nu]);
    Cpmr=permute(Cpmr,[3 2 1]);
    Cps=sqrt(squeeze(sum(sum((C-Cpmr).^2,1),2))./(p.Nr*Nu));
    pred.Cm=Cpm;
    pred.Cs=Cps;
    if nothier
        pred.Cprim=mean(Cpr,2);
        pred.Cpris=std(Cpr,1,2);
    end
end


end