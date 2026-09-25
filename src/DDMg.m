function D=DDMg(t,P,T,a,b,c,p)
%% DDMG: A simple degree day glacier model
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
Dj=zeros(Ne,1);
dTrs=Tr-Ts;

for j=1:Nt
    Tj=T(j)-k2c;
    Tj=Tj-b;
    ddj=Tj-Tm;
    aj=a;
    aj(Dj<0)=a(Dj<0).*(1/0.7); % Based on PyGEM
    Mj=max(aj.*ddj,0);
    frj=(Tj-Ts)./(dTrs);
    frj=max(frj,0);
    frj=min(frj,1);
    fsj=1-frj;
    Sj=(c.*fsj).*(P(j));
    Dj=Dj+Sj-Mj;
    D(j,:)=Dj;   
end


end
