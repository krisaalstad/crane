function r=CRA(p,f,o,hyper,r,g)
%% CRA: Cryospheric ReAnalysis
% Ensemble-based reanalysis with fixed hyperparameters using a simple
% degree day model (DDMg=glacier or DDMs=snow).
% Modify outputs to be EM compatible, need the particle weights and
% parameters for each year (would need residuals for tuning noise term).
% For simplicty just make this the resampled particles which will have
% equal weight. 

% g is the generation counter, just set to one if this is a
% non-hierarchical reanalysis

% NB: The EM part could be accelerated by reusing the proposal (and
% residuals) from the last run, but just adjusting the hyperparameters in 
% AMIS (or similar) to then only iterate if Neff<threshold

% Might be better to output the weighted ensemble and the corresponding
% Gaussian proposal, since this could then be RE-used (for free!) in a
% DM importance sampling scheme for other hyperparameters.

wys=p.wys;
Nwy=numel(wys);
t=f.t; % Synchronize to water years

tstart=datenum(sprintf('%s-%d',p.dstart,p.wys(1)-1));
tend=datenum(sprintf('%s-%d',p.dend,p.wys(end)));
if min(t)>tstart||max(t)<tend
    error('Forcing time series too short')
elseif tstart==t(1)&&tend==t(end)
    % Do nothing
else
    here=t>=tstart&t<=tend;
    t=t(here);
    f.t=f.t(here);
    f.T=f.T(here);
    f.P=f.P(here);
end



Nt=numel(t);
Ne=p.Ne;
Na=p.Na;
mua=hyper.mua; taua=hyper.taua; 
mub=hyper.mub; taub=hyper.taub;
muc=hyper.muc; tauc=hyper.tauc; 
if ~p.glacier
    muv=hyper.muv; tauv=hyper.tauv; 
    p.Np=4;
end
Np=p.Np;

if p.resample
    Nr=p.Nr;
    reh.logZ=0;
    if ~p.samplingchains
        reh.Dpri=zeros(Nt,Nr);
        reh.Dpos=zeros(Nt,Nr);
        reh.Fpri=zeros(Nt,Nr);
        reh.Fpos=nan(Nt,Nr);
        reh.theta=zeros(Nwy,p.Np,Nr);
        reh.theta0=zeros(Nwy,p.Np,Nr);
    end
else
    reh=[];
end

% Check if old runs exist in r for MAGPIES if we want to recycle
if ~isempty(r)&&p.recycle
    % Number of old generations that exist
    Nold=g-1; % Based on generation counter from the "p" struct
else
    Nold=0;
end
logit=@(xb,bd) log(((xb-bd(1))./(bd(2)-bd(1)))...
    ./(1-((xb-bd(1))./(bd(2)-bd(1)))));
expit=@(yu,bd) bd(1)+(bd(2)-bd(1))./(1+exp(-yu));

mu0=[mua; mub; muc]; 
C0=diag([expit(taua,p.sigbd)^2; expit(taub,p.sigbd)^2; expit(tauc,p.sigbd)^2]);
if ~p.glacier
    mu0=[mu0; muv];
    C0=diag([diag(C0); expit(tauv,p.sigbd)^2]);
end
ab=p.ab;
cb=p.cb;
vb=p.vb;

% Force assimF to false if glacier case
if p.glacier
    p.assimF=0;
    AMBinf=2; % Annual mass balance obs error inflation
end

% Array to store the prior proposals used in MAGPIES
theta0=zeros(Np,Ne);

if p.docontrol==1
    a0=5;
    b0=0;
    c0=1;
    if ~p.glacier
        v0=0.5;
    end
    p.Ne=1;
    p.Nre=1;
    p.lowy=1;
    p.lowys=p.wys;
    p.Na=0;
end

sequential=0;
if any(strcmp(p.case,'CP'))||any(strcmp(p.case,'CC'))
    static=1;
else
    static=0;
    if any(strcmp(p.case,'SP'))||any(strcmp(p.case,'AP'))
        sequential=1;
    end
end

logZ=0;

if ~p.samplingchains
    Dpo=zeros(Nt,p.Ne);
    Dpr=zeros(Nt,p.Ne);
    if ~p.glacier
        Fpr=zeros(Nt,p.Ne);
        Fpo=zeros(Nt,p.Ne);
    end
end
theta=zeros(Np,Ne);
thetapo=zeros(numel(wys),Np,Ne); % Output trans posterior params for each year
if ~static
    for k=1:Nwy
        wyk=wys(k);
        tkstart=datenum(sprintf('%s-%d',p.dstart,wyk-1));
        tkend=datenum(sprintf('%s-%d',p.dend,wyk));
        tk=tkstart:tkend;
        tk=tk';
        Ntk=numel(tk);
        herek=t>=min(tk)&t<=max(tk);
        Pk=f.P(herek);
        Tk=f.T(herek);
        hereok=o.t>=min(tk)&o.t<=max(tk);
        to=o.t(hereok);
        Dok=o.D(hereok);
        % Filter out missing values in observations
        missingD=isnan(Dok);
        if any(missingD)
            Dok=Dok(~missingD);
            to=to(~missingD);
        end
        Nok=numel(to);
        if any(to)
            theseo=round(to-tkstart+1); % = Days of WY with obs in tk
            assimDk=1;
        else
            theseo=[];
            assimDk=0;
        end


        assimFk=0;
        if ~p.glacier
            hereosk=o.ts>=min(tk)&o.ts<=max(tk);
            tos=o.ts(hereosk);
            Fosk=o.Fs(hereosk);
            Nosk=numel(tos); 
            missingF=isnan(Fosk);
            if any(missingF)
                Fosk=Fosk(~missingF);
                tos=tos(~missingF);
            end
            if any(tos)
                theseos=round(tos-tkstart+1);
                assimFk=1;
            else
                theseos=[];
                assimFk=0;
            end
            R=[];
            if p.assimD
                R=[R; (p.sigD.^2).*ones(Nok,1)];
            end
            if p.assimF
                R=[R; (p.sigF.^2).*ones(Nosk,1)];
            end
        elseif any(theseo)
            % For glaciers assume annual MB is twice as uncertain
            % as winter MB to account for timing uncertainty.
            R=[p.sigD.^2; (AMBinf.*p.sigD).^2];
        end

        if p.lowy==1&&any(p.lowys==wyk)
            doassimk=0; % Instead of removing obs just turn off DA
        else
            % Check if desired observations actually exist for wyk
            if (p.assimF&&assimFk)||(p.assimD&&assimDk)
                doassimk=1;
            else
                doassimk=0;
            end
        end

        if doassimk==1
            Nak=Na;
            if p.doMAGPIES % This works for both PP and NP
                if Nold>0
                    pros=r(k).pros;
                    prom=r(k).prom;
                    proc=r(k).proc;
                    D=r(k).D;
                    if ~p.glacier
                        F=r(k).F;
                    end
                end
                Nq=(Nold+1)*(Na+1);
                r(k).pros=zeros(Np,Ne,Nq);
                r(k).prom=zeros(Np,Nq);
                r(k).proc=zeros(Np,Np,Nq);
                if ~p.samplingchains
                    r(k).D=zeros(Ntk,Ne,Nq);
                end
                Nqold=Nold*(Na+1);
                if Nold>0
                    r(k).pros(:,:,1:Nqold)=pros;
                    r(k).prom(:,1:Nqold)=prom;
                    r(k).proc(:,:,1:Nqold)=proc;
                    if ~p.samplingchains
                        r(k).D=zeros(Ntk,Ne,Nq);
                    end
                end
                % Allocate prior as first proposal this round
                r(k).prom(:,Nqold+1)=mu0; 
                r(k).proc(:,:,Nqold+1)=C0; 
               
                if ~p.glacier
                    if ~p.samplingchains
                        r(k).F=zeros(Ntk,Ne,Nq);
                    end
                    if Nold>0
                        if ~p.samplingchains
                            r(k).F(:,:,1:Nqold)=F;
                        end
                    end
                end
            end
        else
            Nak=0;
            Nqold=(Nold+1)*(Nak+1); % Not assimilation year, only save prior.
            r(k).D=zeros(Ntk,Ne);
            if ~p.glacier
                r(k).F=zeros(Ntk,Ne);
            end
        end

        for ell=0:Nak
            % Generate prior
            if ell==0&&(sequential==0||k==1)
                alpha0=mua+expit(taua,p.sigbd).*randn(1,Ne);
                a=expit(alpha0,ab);
                beta0=mub+expit(taub,p.sigbd).*randn(1,Ne);
                b=beta0;
                gamma0=muc+expit(tauc,p.sigbd).*randn(1,Ne);
                c=expit(gamma0,cb);
                if ~p.glacier
                    nu0=muv+expit(tauv,p.sigbd).*randn(1,Ne);
                    v=expit(nu0,vb);
                end

                % Assign prior proposal to r for MAGPIES
                if p.doMAGPIES
                    theta0(1,:)=alpha0;
                    theta0(2,:)=beta0;
                    theta0(3,:)=gamma0;
                    if ~p.glacier
                        theta0(4,:)=nu0;
                    end
                    r(k).pros(:,:,Nqold+1)=theta0;
                    r(k).wy=wyk;
                    theta=theta0;
                end
                    

    
             
                
                %
            elseif ell==0&&sequential==1
                if any(strcmp(p.case,'AP'))
                    error('Work in progress :)');
                end
            end
            if p.glacier
                Dk=DDMg(tk,Pk,Tk,a',b',c',p);
            else
                [Dk,Fk]=DDMs(tk,Pk,Tk,a',b',c',v',p);
            end

            if p.doMAGPIES
                if doassimk
                    r(k).D(:,:,Nqold+(ell+1))=Dk;
                    if ~p.glacier
                        r(k).F(:,:,Nqold+(ell+1))=Fk;
                    end
                else
                    r(k).D=Dk;
                    if ~p.glacier
                        r(k).F=Fk;
                    end
                end
                r(k).yk=[]; % Set obs to empty in case no assimilation year
            end

            if ell==0
                Dpr(herek,:)=Dk;
                if ~p.glacier
                    Fpr(herek,:)=Fk;
                end
            end

            if ell<Na&&doassimk
                Dpk=Dk(theseo,:);
                % Clamping to avoid precision issues at boundaries
                epsis=1e-3;
                a(a<=(ab(1)+epsis))=ab(1)+epsis;
                a(a>=(ab(2)-epsis))=ab(2)-epsis;
                theta(1,:)=logit(a,ab);
                theta(2,:)=b;
                c(c<=(cb(1)+epsis))=cb(1)+epsis;
                c(c>=(cb(2)-epsis))=cb(2)-epsis;
                theta(3,:)=logit(c,cb);
                if p.glacier
                    if ~p.bsa
                        theta=EnKA(theta,Dok,Dpk,Na,R,p.dostoch);
                    else % Assimilate summer and not net MB.
                        % ba=bw-bs <-> bs=bw-ba
                        Doka=Dok;
                        Doka(2)=Doka(1)-Doka(2);
                        Dpka=Dpk;
                        Dpka(2)=Dpka(1)-Dpka(2);
                        theta=EnKA(theta,Doka,Dpka,Na,R,p.dostoch);
                    end
                else
                    v(v<=(vb(1)+epsis))=vb(1)+epsis;
                    v(v>=(vb(2)-epsis))=vb(2)-epsis;
                    theta(4,:)=logit(v,vb);
                    Fpk=Fk(theseos,:);
                    if p.assimF&&p.assimD
                        Ypk=[Dpk; Fpk];
                        yk=[Dok; Fosk];
                    elseif p.assimD&&~p.assimF
                        Ypk=Dpk;
                        yk=Dok;
                    elseif p.assimF&&~p.assimD
                        Ypk=Fpk;
                        yk=Fosk;
                    end
                    theta=EnKA(theta,yk,Ypk,Na,R,p.dostoch);
                end
                % Proposal construction
                if (ell==(Na-1))||(p.doMAGPIES==1)
                    propm=mean(theta,2);
                    anom=theta-propm;
                    propc=(anom*anom')./Ne;
                    try
                        L=chol(propc,'lower');
                    catch
                        propc=NSPD(propc);
                        L=chol(propc,'lower');
                    end
                    z=randn(Np,Ne);
                    theta=propm+L*z;
                    prop=theta;
                    if p.doMAGPIES
                        qis=Nqold+ell+2; % Iteration counter in MAGPIES
                        r(k).pros(:,:,qis)=prop;
                        r(k).prom(:,qis)=propm;
                        r(k).proc(:,:,qis)=propc;
                    end
                end
                a=expit(theta(1,:),ab);
                b=theta(2,:);
                c=expit(theta(3,:),cb);
                if ~p.glacier
                    v=expit(theta(4,:),vb);
                end
            elseif ell==p.Na&&doassimk==1 % PIES or MAGPIES update
                if p.doMAGPIES
                    Dpk=r(k).D(theseo,:,:);
                else
                    Dpk=Dk(theseo,:);
                end
                if ~p.glacier
                    if p.doMAGPIES
                        Fpk=r(k).F(theseos,:,:);
                    else
                        Fpk=Fk(theseos,:);
                    end
                    if p.assimF&&p.assimD
                        Ypk=[Dpk; Fpk];
                        yk=[Dok; Fosk];
                        r(k).theseD=theseo;
                        r(k).theseF=theseos;
                    elseif p.assimD&&~p.assimF
                        Ypk=Dpk;
                        yk=Dok;
                        r(k).theseD=theseo;
                        r(k).theseF=[];
                    elseif p.assimF&&~p.assimD
                        Ypk=Fpk;
                        yk=Fosk;
                        r(k).theseF=theseos;
                        r(k).theseD=[];
                    end
                else
                    Ypk=Dpk;
                    yk=Dok;
                    if p.bsa
                        % Assim summer mass balance not net
                        % ba=bw-bs <-> bs=bw-ba
                        yk(2)=yk(1)-yk(2);
                        Ypk(2,:,:)=Ypk(1,:,:)-Ypk(2,:,:);
                    end
                    if p.doMAGPIES
                        r(k).theseD=theseo;
                    end
                end
                
                % Obs are always oredered D first then F (if both exist).
                % Recall that MAGPIES can be computed for new psi outside
                % this routine based on the outputs saved in "r"
                if p.doMAGPIES
                    %error('MAGPIES')
                    [w,logZk]=MAGPIES(yk,Ypk, R,...
                        mu0, diag(C0), r(k).pros, r(k).prom, r(k).proc );
                    r(k).w=w;
                    r(k).logZ=logZk;
                    Neff=1./sum(w(:).^2);
                    r(k).Neff=Neff;
                    r(k).yk=yk; % Save obs if final iteration and assimilation year.
                    

                else % PIES
                    [w,logZk] = PIES(yk,Ypk,R,...
                        mu0,C0,prop,propm,propc); 
                    logZ=logZ+logZk;
                    Neff=1./sum(w.^2);
                    reinds=randsample(Ne,p.Nre,true,w);
                    Dk=Dk(:,reinds);
                    a=a(reinds);
                    b=b(reinds);
                    c=c(reinds);
                    Dpo(herek,:)=Dk;
                    if ~p.glacier
                        Fk=Fk(:,reinds);
                        v=v(reinds);
                        Fpo(herek,:)=Fk;
                    end
                    thetapo(k,:,:)=theta(:,reinds);
                end
            elseif doassimk==0
                Dpo(herek,:)=Dk;
                if ~p.glacier
                    Fpo(herek,:)=Fk;
                end
                logZk=0;
                logZ=logZ+logZk;
            end


            if p.resample&&ell==Nak
                if p.Ng>1
                    error('Resampled history currently only for Ng=1');
                end
                here=f.t>=tkstart&f.t<=tkend;
                herewy=wys==wyk;

                % Output a resampled history
                prire=randsample(p.Ne,Nr,true); % Resample prior
                if Nr==p.Ne
                    % Special case, no resampling needed for prior
                    prire=1:p.Ne;
                end
                Dk=r(k).D;
                reh.Dpri(here,:)=Dk(:,prire,1); 
                if ~p.glacier
                    Fk=r(k).F;
                    reh.Fpri(here,:)=Fk(:,prire,1);
                end
                
                if doassimk==1
                    reh.logZ=reh.logZ+logZk;
                    Nh=p.Ne*(p.Na+1);
                    posre=randsample(Nh,Nr,true,w(:));
                    if ~p.samplingchains
                        Dk=reshape(Dk,size(Dk,1),Nh);
                        reh.Dpos(here,:)=Dk(:,posre);
                        if ~p.glacier
                            Fk=reshape(Fk,size(Fk,1),Nh);
                            reh.Fpos(here,:)=Fk(:,posre);
                        end
                        thetah=reshape(r(k).pros,p.Np,Nh);
                        reh.theta(herewy,:,:)=thetah(:,posre);
                        reh.theta0(herewy,:,:)=theta0(:,prire);
                    end
                else
                    if ~p.samplingchains
                        reh.Dpos(here,:)=reh.Dpri(here,:);
                        if ~p.glacier
                            reh.Fpos(here,:)=reh.Fpri(here,:);
                        end
                        reh.theta(herewy,:,:)=theta0(:,prire);
                        reh.theta0(herewy,:,:)=theta0(:,prire);
                    end
                end
            end
                        
          
        end
    end

else % Complete pooling or climatological (static global parameters) 
    fprintf('\n Static case=%s \n',p.case);

    
    if p.doMAGPIES
        Dh=zeros(Nt,p.Ne,p.Na+1);
        if ~p.glacier
            Fh=zeros(Nt,p.Ne,p.Na+1);
        end
    end

    for ell=0:Na
        Do=[];
        Dp=[];
        if ~p.glacier
            Fo=[];
            Fp=[];
        end
        if p.doMAGPIES&&(ell==0)
            prop=zeros(p.Np,p.Ne,p.Na+1);
            propm=zeros(p.Np,p.Na+1);
            propc=zeros(p.Np,p.Np,p.Na+1);
            propm(:,1)=mu0;
            propc(:,:,1)=C0;
        end
        % Generate prior
        if ell==0&&(sequential==0)%||k==1
            a=expit(mua+expit(taua,p.sigbd).*randn(1,Ne),ab);
            b=mub+expit(taub,p.sigbd).*randn(1,Ne);
            c=expit(muc+expit(tauc,p.sigbd).*randn(1,Ne),cb);
            if ~p.glacier
                v=expit(muv+expit(tauv,p.sigbd).*randn(1,Ne),vb);
            end
            if p.docontrol
                a=a0;
                b=b0;
                c=c0;
                v=v0;
            end

            if p.doMAGPIES&&ell==0 % ell=0 just for clarity
                prop(1,:,ell+1)=logit(a,ab);
                prop(2,:,ell+1)=b;
                prop(3,:,ell+1)=logit(c,cb);
                if ~p.glacier
                    prop(4,:,ell+1)=logit(v,vb);
                end
                theta0=squeeze(prop(:,:,ell+1)); % Save prior parameters
            end
            %
        end
        for k=1:numel(wys)
            wyk=wys(k);
            tkstart=datenum(sprintf('%s-%d',p.dstart,wyk-1));
            tkend=datenum(sprintf('%s-%d',p.dend,wyk));
            tk=tkstart:tkend;
            tk=tk';
            Ntk=numel(tk);
            if (p.lowy==1&&any(p.lowys==wyk))||tk(1)>max(o.t)
                doassimk=0; % Instead of removing obs just turn off DA
            else
                doassimk=1; 
            end
            herek=t>=min(tk)&t<=max(tk);
            Pk=f.P(herek);
            Tk=f.T(herek);
            hereok=o.t>=min(tk)&o.t<=max(tk);
            to=o.t(hereok);
            Dok=o.D(hereok);
            Nok=numel(to);
            theseo=round(to-tkstart+1); % = Days of WY with obs in tk
            if ~p.glacier
                hereosk=o.ts>=min(tk)&o.ts<=max(tk);
                tos=o.ts(hereosk);
                Fosk=o.Fs(hereosk);
                Nosk=numel(tos);
                theseos=round(tos-tkstart+1);
            end


            if p.glacier
                Dk=DDMg(tk,Pk,Tk,a',b',c',p);
            else
                [Dk,Fk]=DDMs(tk,Pk,Tk,a',b',c',v',p);
            end
            Dpk=Dk(theseo,:);
            if ~p.glacier
                Fpk=Fk(theseos,:);
            end
            if doassimk==1
                Do=[Do; Dok];
                Dp=[Dp; Dpk];
                if ~p.glacier
                    Fo=[Fo; Fosk];
                    Fp=[Fp; Fpk];
                end
            end
            if ell<Na
                if ell==0
                    Dpr(herek,:)=Dk;
                    if ~p.glacier
                        Fpr(herek,:)=Fk;
                    end
                end
            elseif ell==Na
                Dpo(herek,:)=Dk;
                if ~p.glacier
                    Fpo(herek,:)=Fk;
                end
            end
            if p.doMAGPIES
                Dh(herek,:,ell+1)=Dk;
                if ~p.glacier
                    Fh(herek,:,ell+1)=Fk;
                end
            end
        end
    
        if p.resample&&(ell==0)
            prire=randsample(p.Ne,p.Nr,true); % Resample prior
            if Nr==p.Ne
                % Special case, no resampling needed for prior
                prire=1:p.Ne;
            end
            reh.Dpri=squeeze(Dh(:,prire,1));
            if ~p.glacier
                reh.Fpri=squeeze(Fh(:,prire,1)); % Changed from Fpr(:,prire,1)
            end
        end


        if ell<Na&&~p.docontrol
            theta(1,:)=logit(a,ab);
            theta(2,:)=b;
            theta(3,:)=logit(c,cb);
          
            if p.glacier
                Yp=Dp;
                y=Do;
                if p.bsa
                    % ba=bw-bs <-> bs=bw-ba
                    y(2:2:end)=y(1:2:(end-1))-y(2:2:end);
                    Yp(2:2:end,:)=Yp(1:2:(end-1),:)-Yp(2:2:end,:);
                end

                if ell==0
                    R=ones(numel(Do),1);
                    R(1:2:(end-1))=p.sigD.^2;
                    R(2:2:end)=(AMBinf.*p.sigD).^2;
                end
            else
                theta(4,:)=logit(v,vb);
                if ell==0
                    R=[];
                end
                if p.assimF&&p.assimD
                    Yp=[Dp; Fp];
                    y=[Do; Fo];
                    if ell==0
                        R=[R; (p.sigD.^2).*ones(numel(Do),1);...
                            (p.sigF.^2).*ones(numel(Fo),1)];
                    end
                elseif p.assimD&&~p.assimF
                    Yp=Dp;
                    y=Do;
                    if ell==0
                        R=[R; (p.sigD.^2).*ones(numel(Do),1)];
                    end
                elseif p.assimF&&~p.assimD
                    Yp=Fp;
                    y=Fo;
                    if ell==0
                        R=[R; (p.sigF.^2).*ones(numel(Fo),1)];
                    end
                end
                if numel(R)~=numel(y)
                    error('stop!');
                end
            end
            theta=EnKA(theta,y,Yp,Na,R,0);
            if ell==(Na-1)&&~(p.doMAGPIES)
                propm=mean(theta,2);
                anom=theta-propm;
                propc=(anom*anom')./Ne;
                L=chol(propc,'lower');
                z=randn(Np,Ne);
                theta=propm+L*z;
                prop=theta;
            elseif p.doMAGPIES
                ellis=ell+2;
                propm(:,ellis)=mean(theta,2);
                anom=theta-propm(:,ellis);
                propc(:,:,ellis)=(anom*anom')./Ne;
                L=chol(propc(:,:,ellis),'lower');
                z=randn(Np,Ne);
                theta=propm(:,ellis)+L*z;
                prop(:,:,ellis)=theta;

                % Create prediction array across proposals
                if ell==0
                    Ypp=zeros(size(Yp,1),size(Yp,2),Na+1);
                end
                % Stores predictions (from this iteration, not yet updated)
                Ypp(:,:,ell+1)=Yp;
            end
            a=expit(theta(1,:),ab);
            b=theta(2,:);
            c=expit(theta(3,:),cb);
            if ~p.glacier
                v=expit(theta(4,:),vb);
            end
        elseif ~p.docontrol % PIES update
            if ~p.glacier
                if p.assimF&&p.assimD
                    Yp=[Dp; Fp];
                    y=[Do; Fo];
                elseif p.assimD&&~p.assimF
                    Yp=Dp;
                    y=Do;
                elseif p.assimF&&~p.assimD
                    Yp=Fp;
                    y=Fo;
                end
            else
                Yp=Dp;
                y=Do;
            end
            if ~p.doMAGPIES
                w = PIES(y,Yp,R,...
                    mu0,C0,prop,propm,propc);
                Neff=1./sum(w.^2);
                reinds=randsample(Ne,Ne,true,w);
                Dk=Dk(:,reinds);
                a=a(reinds);
                b=b(reinds);
                c=c(reinds);
                Dpo=Dpo(:,reinds);
                if ~p.glacier
                    Fpo=Fpo(:,reinds);
                    v=v(reinds);
                end
            else
                if Na>0
                    Ypp(:,:,end)=Yp; % Add latest prediction from ell=Na
                    [w,~]=MAGPIES(y,Ypp,R,mu0,diag(C0),...
                        prop,propm,propc);
                    reh.Neff=1./sum(w(:).^2);
                    

                    % Resample history based on weights
                    % Resample should be activated always with CP and NP
                    if p.resample
                        Nh=p.Ne*(p.Na+1);
                        posre=randsample(Nh,Nr,true,w(:));
                        Dh=reshape(Dh,Nt,Nh);
                        reh.Dpos=Dh(:,posre);
                        if ~p.glacier
                            reh.Fpos=Fh(:,posre);
                        end
                        thetah=reshape(prop,p.Np,Nh);
                        reh.theta=thetah(:,posre);
                        reh.theta0=theta0;
                    end
                end
            end
        end
    end
end

% Now just outputs a single struct r containing the results
if ~p.doMAGPIES
    % Dpo,Dpr,logZ, thetapo
    r.Dpo=Dpo;
    r.Dpr=Dpr;
    r.logZ=logZ;
    r.thetapo=thetapo;
elseif p.resample
    % Only output the resampled history if this was flagged.
    r=reh;
    r.t=t;
    clear reh;
end
% For MAGPIES this structure has already been created and might contain
% several generations of results.





end


%{



tic;
Dclast=0;
for k=1:numel(wys)
    wyk=wys(k);
    tk=datenum(sprintf('01-Sep-%d',wyk-1)):...
        datenum(sprintf('31-Aug-%d',wyk));
    tk=tk';
    herek=t>=min(tk)&t<=max(tk);
    Pk=f.P(herek);
    Tk=f.T(herek);
    Dk=DDM(tk,Pk,Tk,a,b,c);
    %{
    if wyk==1951D
        Dck=Dk+Dlast;
    else
        Dck=Dk;
    end
    %}
    D(herek,:)=Dk;
    Dck=Dk+Dclast;
    Dc(herek,:)=Dck;
    Dclast=Dck(end,:); % CSMB at end of the year
    D(herek,:)=Dk;

    %{
    figure(1);
    subplot(1,3,1);
    plot(tk,cumsum(max(3.*(Tk-273.15),0)));
    subplot(1,3,2);
    plot(tk,cumsum(Pk));
    subplot(1,3,3);
    plot(tk,Dk);
    datetick;
    return
    %}
end
toc;
%}
