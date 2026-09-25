function [D,fSCA]=DDMs(t,P,T,a,b,c,v,p)
%% DDMs: A simple degree day seasonal snow model
% Runs potentially large ensemble for a single point/cell

% Hyperparameters:
Tr=p.Tr; % Temperature threshold for rainfall (degrees C)
Ts=p.Ts; % Temperature threshold for snowfall (degrees C)
Tm=p.Tm; % Temperature threshold for snowmelt (degrees C)
k2c=273.15;


% Pre-allocation
Nt=numel(t); % Number of time steps.
Ne=numel(a);
D=zeros(Nt,Ne);
fSCA=zeros(Nt,Ne);
Dj=zeros(Ne,1);
fSCAj=zeros(Ne,1);
dTrs=Tr-Ts;


[mu,Dm,nada]=deal(zeros(Ne,1));
issnow=mu>0;

for j=1:Nt
    Tj=T(j)-k2c;
    Tj=Tj-b;
    ddj=Tj-Tm;
    Mj=max(a.*ddj,0);
    frj=(Tj-Ts)./(dTrs);
    frj=max(frj,0);
    frj=min(frj,1);
    fsj=1-frj;
    Sj=(c.*fsj).*(P(j));
    NAj=Sj-Mj;


    % Track current melt depth.
    Dmcur=Dm;

    % Update accumulated melt SWE depth
    Dm=max(Dm-NAj,0).*issnow; 

    % Update peak SWE depth
    mu=mu+max(NAj-Dmcur,0); % Accounts for residual melt depth

    % Flag if a snowpack has been created
    issnow=mu>0;

    if any(issnow)
        [fSCAj(issnow),Dj(issnow)]=fSCA_scheme(mu(issnow),Dm(issnow),v(issnow));
        fSCAj(~issnow)=nada(~issnow);
        Dj(~issnow)=nada(~issnow);
    else
        fSCAj=nada;
        Dj=nada;
    end
    
    % Reset accounting variables if snowpack is gone
    issnow=fSCAj>0;
    Dm(~issnow)=0;
    mu(~issnow)=0;

    D(j,:)=Dj;
    fSCA(j,:)=fSCAj;
end




    %% fSCA scheme.
    function [Fs,Das]=fSCA_scheme(mus,Dms,chis)
        % Based on Liston (2004; JClim).
        % Calculates the fSCA given an accumulated melt depth increment 
        % D_m_new-D_m_cur and an initial (premelt) lognormal snow distribution
        % with mean mu_new and coefficient of variation (CV=standard dev/mean).
            sdt=sqrt(log(1+chis.^2));
            mut=log(mus)-0.5.*sdt.^2;
            mut(mut==-Inf)=0;
            z=(log(Dms)-mut)./(sdt);

            Fs=0.5.*erfc(z./sqrt(2));
            Fs(Fs<1e-2)=0;
            Das=0.5.*exp(mut+0.5.*sdt.^2).*...
                erfc((z-sdt)./sqrt(2))...
                -Fs.*Dms;
            Das(Das<0)=0;
            Das(Fs==0)=0;
            
    end

end
