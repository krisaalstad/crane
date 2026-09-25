function gen = particleEM(p,f,o)
% Wrapper function for running multiple generations of the particle EM
% routine based on particle-based inference via MAGPIES.

% Gradients and EM objective f
fevalp=@(w,thens,mu,sig,tau,m,s,eta,chi) -(0.5.*((mu-m)./s).^2)...
    -(0.5.*((tau-eta)./chi).^2)...
    -sum(w.*((((thens-mu).^2)./(2.*sig.^2))+log(sig)));
dfdmu=@(w,thens,mu,sig,m,s) ((m-mu)./(s^2))...
    +(sum(w.*(thens-mu))./(sig.^2));
d2fdmu2=@(w,thens,sig,s) (-1/(s^2))...
    -(sum(w.*ones(size(thens)))./(sig.^2));
dsigdtau=@(sig,bd) (sig-bd(1)).*(bd(2)-sig)./(bd(2)-bd(1));
d2sigdtau2=@(sig,dsdt,bd) (bd(2)+bd(1)-2.*sig).*dsdt./(bd(2)-bd(1));
dfdtau=@(w,thens,mu,sig,dsdt,tau,eta,chi) ((eta-tau)./(chi^2))...
    +sum(w.*(((((thens-mu).^2)./(sig.^3))-(1./sig)).*dsdt));
d2fdtau2=@(w,thens,mu,sig,dsdt,d2sdt2,chi) (-1/(chi^2))...
    +(sum(w.*(...
    (((1./sig.^2)-3.*((thens-mu).^2)./(sig.^4)).*(dsdt).^2)+...
    ((((thens-mu).^2)./(sig.^3))-(1./sig)).*d2sdt2)));
d2fdtaudmu=@(w,thens,mu,sig,dsdt) -2.*(sum(w.*(((thens-mu)./(sig.^3)).*dsdt)));

 

Ng=p.Ng;
maxit=p.EMit;
psi=zeros(2*p.Np,maxit+1); % psi vector for updating with EM
gen.r=[]; % Initialize with empty result struct
gen.psi=zeros(p.Np*2,maxit+1,Ng);
gen.p=p;
Nwy=numel(p.wys);

fval=zeros(maxit+1,Ng);



for g=1:Ng
    gen.g=g; % Update generation counter in "gen" struct to track generations
    if g==1
        flast=-Inf; % Should ideally be re-sampled (new seed) at each new M-step to account for sampling approximation.
        % But we don't do this at least when recycling, since that defeats
        % the point of recycling.
    end
    for k=1:maxit
        % Start at the hyperprior mean, possibly perturbed.
        if k==1
            if g==1
                pinit=0; % Perturbed initialization, true or false
                h.mua=p.ma+pinit.*p.sa.*randn(1,1);
                h.mub=p.mb+pinit.*p.sb.*randn(1,1);
                h.muc=p.mc+pinit.*p.sc.*randn(1,1);

                h.taua=p.ea+pinit.*p.xa.*randn(1,1);
                h.taub=p.eb+pinit.*p.xb.*randn(1,1);
                h.tauc=p.ec+pinit.*p.xc.*randn(1,1);
                if ~p.glacier
                    h.muv=p.mv+pinit.*p.sv.*randn(1,1);
                    h.tauv=p.ev+pinit.*p.xv.*randn(1,1);
                end


                psi(1,1)=h.mua;
                psi(2,1)=h.mub;
                psi(3,1)=h.muc;
                psi(p.Np+1,1)=h.taua;
                psi(p.Np+2,1)=h.taub;
                psi(p.Np+3,1)=h.tauc;
                if ~p.glacier
                    psi(p.Np,1)=h.muv;
                    psi(2*p.Np,1)=h.tauv;
                end


                % m: Hyperprior mean on prior mean mu
                hm=p.m;

                % s: Hyperprior std on prior mean mu
                hs=p.s;

                % eta: Hyperprior mean on prior scale tau
                he=p.e;

                % chi: Hyperprior mean on prior scale tau
                hx=p.x;
            else
                % Use best estimate from last generation if start of new
                % generation
                psi(:,:)=0; % For hygiene, not strictly needed.
                
                psi(:,1)=gen.psi(:,end,g-1); % Best estimate form last generation
               
            end

            
        end

        if k==1||~p.recycle
             % Use updated hyperparameters from last generation in CRA
             h.mua=psi(1,k);
             h.mub=psi(2,k);
             h.muc=psi(3,k);
             h.taua=psi(p.Np+1,k);
             h.taub=psi(p.Np+2,k);
             h.tauc=psi(p.Np+3,k);
             if ~p.glacier
                 h.muv=psi(p.Np,k);
                 h.tauv=psi(2*p.Np,k);
             end


            % Run reanalysis only for the first EM k=1 iteration of a given
            % generation (g). Or always if you don't want to recycle.
            gen.r=CRA(p,f,o,h,gen.r,gen.g);
            obsyear=zeros(Nwy,1,'logical');
            if p.recycle
                Nq=g*(p.Na+1);
            else
                Nq=p.Na+1;
            end
            Nh=p.Ne*Nq;
            if ~p.resample
                theta=zeros(p.Np,p.Ne,Nq,Nwy); % Concatenate later to "particle history"
                w=zeros(p.Ne,Nq,Nwy);
            else
                theta=zeros(p.Np,p.Ne,Nwy); % Concatenate later to resampled "particle history"
                w=zeros(p.Ne,Nwy);
                weq=(1/p.Ne).*ones(p.Ne,1);
            end
            for yy=1:Nwy
                if ~isempty(gen.r(yy).yk)
                    wn=gen.r(yy).w;
                    thetan=gen.r(yy).pros;
                    obsyear(yy)=1;
                    if p.resample
                        resample=randsample(Nh,p.Ne,true,wn(:));
                        thetan=reshape(thetan,p.Np,Nh);
                        thetan=thetan(:,resample);
                        theta(:,:,yy)=thetan;
                        wn=weq;
                        w(:,yy)=wn;
                    else
                        w(:,:,yy)=wn;
                        theta(:,:,:,yy)=thetan;
                    end
                end
            end
            Nwyo=sum(obsyear); % Number of wys with obs
            if p.resample
                w=w(:,obsyear);
                theta=theta(:,:,obsyear);
            else
                w=w(:,:,obsyear); % <- Can do new updates for w for new psi
                theta=theta(:,:,:,obsyear);
                theta=reshape(theta,[p.Np,Nh,Nwyo]);
            end
            wc=w(:); % Concatenate weights a vector
        else
            % Recycle old particles:
            % Do reMAGPIES for later EM iterations
            [w,~,thetare] = reMAGPIES(p,psi(:,k),gen);
            if p.resample
                theta=thetare;
            end
            wc=w(:); % Concatenate weights to (Nh x Nwyo) vector
        end
    
        
        lambda=1;
        lambdaf=10;
        if k<=maxit
            optits=20; % Maximization iterations
            psiell=psi(:,k);
            % Note that ell and ll are not the same here (ll is for the lambda
            % adaptation)
            for ell=1:optits
                Jj=zeros(2*p.Np,1);
                Hj=zeros(2*p.Np,2*p.Np);
                fell=0;
                for j=1:p.Np % Loop over post parameter vector entries
                    % Populate Jacobian for each hyperparam associated with each
                    % post param

                    % Equal (resampled) weights for now
                    thetaj=reshape(squeeze(theta(j,:,:)),[],1); % (Nh x Nwy) vector

                    
                    
                    tauj=psiell(p.Np+j);
                    sigj=p.gexpit(tauj,p.sigbd);
                    dsdt=dsigdtau(sigj,p.sigbd);

                    Jj(j)=dfdmu(wc,thetaj,psiell(j),sigj,hm(j),hs(j));
                    Hj(j,j)=d2fdmu2(wc,thetaj,sigj,hs(j));

                    d2sdt2=d2sigdtau2(sigj,dsdt,p.sigbd);
                    Jj(p.Np+j)=dfdtau(wc,thetaj,psiell(j),sigj,dsdt,tauj,he(j),hx(j));
                    Hj(p.Np+j,p.Np+j)=d2fdtau2(wc,thetaj,psiell(j),sigj,dsdt,d2sdt2,hx(j));
                    Hj(j,p.Np+j)=d2fdtaudmu(wc,thetaj,psiell(j),sigj,dsdt);
                    Hj(p.Np+j,j)=Hj(j,p.Np+j); % Hessian is symmetric


                    fell=fell+fevalp(wc,thetaj,psiell(j),sigj,tauj,...
                        hm(j),hs(j),he(j),hx(j));
                end
                if k==1&&ell==1
                    fval(1,g)=fell;
                    flast=-Inf; 
                end
                
                % Old: ~p.recycle&&ell==1&&k==1
                % Now want to reset for each generation
                if (~p.recycle||k==1)&&ell==1
                    % Reset flast to account for Monte Carlo sampling error
                    flast=fell; 
                end

                Jj=-Jj;
                Hj=-Hj;

                Hj=(Hj+Hj')./2;
                try
                    tmp=chol(Hj);
                catch
                    % Saddle free method
                    [V,L]=eig(Hj,'vector');
                    Hj=V*diag(max(abs(L),1e-6))*V';
                end
                Dj=diag(max(diag(Hj),1e-6));
                Ij=Dj;
                psinew=psiell-((Hj+lambda.*Ij)\Jj);
                ll=0;
                llmax=20;
                lambdamax=1e7;
                lambdamin=1e-7;
                while ll<llmax
                    fell=0;
                    for j=1:p.Np
                        thetaj=reshape(squeeze(theta(j,:,:)),[],1); % (Nh x Nwy) vector
                        taunew=psinew(p.Np+j);
                        signew=p.gexpit(taunew,p.sigbd);
                        fell=fell+fevalp(wc,thetaj,psinew(j),signew,taunew,...
                            hm(j),hs(j),he(j),hx(j));

                    end
                    if fell>=flast
                        lambda=lambda/(2.*lambdaf); % Reduce lambda for next ell
                        lambda=max(lambda,lambdamin);
                        ll=llmax;
                        flast=fell;
                    else
                        lambda=lambda*lambdaf; % Increase lambda and try again
                        if lambda>lambdamax
                            lambda=lambdamax;
                            ll=llmax;
                        else
                            ll=ll+1;
                        end
                        psinew=psiell-((Hj+lambda.*Ij)\Jj);
                        if ll==llmax
                            % Failed to improve, reject move and set
                            % psinew to psiell to avoid moving slowly
                            % in the wrong direction.
                            psinew=psiell;
                        end
                    end
                end
                psiell=psinew;
            end
            fval(k+1,g)=flast;
            psi(:,k+1)=psiell;

            if k==maxit&&g==p.Ng % Evaluate Hessian at MAP-II estimate
                Hhat=zeros(2*p.Np,2*p.Np);

                for j=1:p.Np % Loop over post parameter vector entries
                    % Populate Jacobian for each hyperparam associated with each
                    % post param
                    thetaj=reshape(squeeze(theta(j,:,:)),[],1); % (Nh x Nwy) vector

                    tauj=psiell(p.Np+j);
                    sigj=p.gexpit(tauj,p.sigbd);
                    dsdt=dsigdtau(sigj,p.sigbd);
                    d2sdt2=d2sigdtau2(sigj,dsdt,p.sigbd);

                    Hhat(j,j)=d2fdmu2(wc,thetaj,sigj,hs(j));
                    Hhat(p.Np+j,p.Np+j)=d2fdtau2(wc,thetaj,psiell(j),sigj,dsdt,d2sdt2,hx(j));
                    Hhat(j,p.Np+j)=d2fdtaudmu(wc,thetaj,psiell(j),sigj,dsdt);
                    Hhat(p.Np+j,j)=Hhat(j,p.Np+j); % Hessian is symmetric
                end
                Hhat=-Hhat; % Hessian for phi not f
            end
        end
    end
    gen.psi(:,:,g)=psi;
end
gen.Hhat=Hhat; % Final approximation of precision matix from Laplace
plot(fval(:))
gen.fval=fval;








end
