function [w,logZ,thetare] = reMAGPIES(p,psi,gen)
% Rerun MAGPIES for the hyperparameters in psi (size: Np x Ns)

Ns=size(psi,2); % Can be one or more.
Nwy=numel(p.wys);
Nq=gen.g*(p.Na+1);
w=zeros(p.Ne,Nq,Nwy,Ns); % Reset w for new MAGPIES evaluation
if p.resample&&Ns==1
    wre=zeros(p.Ne,Nwy); % Resample Ne (Not Nq*Ne) equal weights
    thetare=zeros(p.Np,p.Ne,Nwy); % Resampled theta for Ns=1 & p.recycle=1
else
    thetare=[];
end
logZ=zeros(Nwy,Ns);
prim=psi(1:p.Np,:);
priv=psi((p.Np+1):(2*p.Np),:);
priv=p.gexpit(priv,p.sigbd).^2;
obsyear=zeros(Nwy,1,'logical');
for n=1:Nwy
    if ~isempty(gen.r(n).yk)
        obsyear(n)=1;
        yn=gen.r(n).yk;
        NoD=numel(gen.r(n).theseD); % 0 if empty []
        No=NoD;
        if ~p.glacier % Glacier runs don't have fSCA obs
            NoF=numel(gen.r(n).theseF);
        else
            NoF=0;
        end
        No=No+NoF;
        pred=zeros(No,p.Ne,Nq);
        R=zeros(No,1);
        if NoD>0
            D=gen.r(n).D(gen.r(n).theseD,:,:);
            pred(1:NoD,:,:)=D;
            R(1:NoD)=p.sigD^2;
            if p.glacier
                AMBinf=2;
                R(2:2:NoD) = (AMBinf.*p.sigD).^2;
            end
        end
        if NoF>0
            F=gen.r(n).F(gen.r(n).theseF,:,:);
            pred((NoD+1):No,:,:)=F;
            R((NoD+1):No)=p.sigF^2;
        end
        pros=gen.r(n).pros;
        prom=gen.r(n).prom;
        proc=gen.r(n).proc;
        
        for k=1:Ns
            [wn,logZ(n,k)]=MAGPIES(yn,pred,R,prim(:,k),priv(:,k),pros,prom,proc);
            w(:,:,n,k)=wn;
        end
        % Recall: theta does not change here only w via psi.

        % Optional resampling step
        if p.resample&&Ns==1
            wn=wn(:);
            Nh=p.Ne*Nq;
            resample=randsample(Nh,p.Ne,true,wn);
            thetaren=reshape(pros,p.Np,Nh); 
            thetaren=thetaren(:,resample);
            thetare(:,:,n)=thetaren;
            wre(:,n)=1/p.Ne;
            % We could also resample to smaller arrays (e.g. Ne) to speed
            % up
        end
    end
end
w=w(:,:,obsyear,:); % <- Can do new updates for w for new psi
if p.resample&&Ns==1
    w=wre(:,obsyear);
    thetare=thetare(:,:,obsyear);
end
logZ=logZ(obsyear,:);

end

