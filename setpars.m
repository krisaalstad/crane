%% Set up run parameters for crane

if any(exist('p','var'))
    clear p;
end
try 
    p.glacier=glacier; % Glacier flag, 1=glacier (allows D<0), 0=seasonal snow
catch
    p.glacier=0;
end

% Select the case:
% Options:
% CP=Complete pooling
% PP=Partial pooling
% NP=No pooling
p.case='PP';

% Select the case:
p.Ng=2; % Number of generations to run the particle EM
p.EMit=10; % Number of EM iterations for each generation
p.iota=p.Ng*p.EMit;
p.recycle=1; % Recycle particles in EM iterations?
p.resample=0; % Resample inner loop for posterior predictions
p.lowy=0;
p.lowys=[];%(2005:2024)';
p.bsa=0; % Assimilate summer balance (1=summer, 0=net annual)
p.samplingchains=0; % Only active for PMCMC sampling, should otherwise be 0.
p.Ncores=20;%21;
p.dopar=0;
p.docontrol=0;
p.dolocal=0;
p.onlya=1; % Climatological Calibration DA
%p.schemes={'Prior','PIES'};
p.Np=3;
if ~p.glacier
    p.Np=p.Np+1;
end
p.Na=4; % Na seems to be the issue (>4), but not sure why?
% Maybe there is a subtle issue with the MAGPIES update?
p.Ne=1e2;
p.doMAGPIES=1;
p.dostoch=0;
p.Nr=1e2; % Number of particles to resample in final inner loop
p.Nss=1e3; % Number of psi subsamples to pick from PMCMC chains



p.dstart='01-Sep';
p.dend='31-Aug';

%% Hyperparameters
p.Tr=3;
p.Ts=-1;
p.Tm=0;
p.ab=[1 7]; 
p.cb=[0.1 5]; 
p.vb=[1e-2 1]; %1e-2 1
p.SCD_type=1;

%% Transforms
p.glogit=@(xb,bd) log(((xb-bd(1))./(bd(2)-bd(1)))...
    ./(1-((xb-bd(1))./(bd(2)-bd(1)))));
p.gexpit=@(yu,bd) bd(1)+(bd(2)-bd(1))./(1+exp(-yu));


% Hyperprior bounds on location hyperparameters are the same as bounds
% on the underlying parameter. 
% Following the notation in Table 2, but using a,b,c,v instead of alpha,
% beta, gamma, nu for economy. So e.g. ma is the same as m_alpha in Table
% 2, mc is the same as m_gamma etc...

% Hyperprior location and scale for the prior location hyperparameters (mu)
p.ma=p.glogit(4,p.ab); % 5
p.mb=0;
p.mc=p.glogit(1,p.cb);
p.mv=p.glogit(0.4,p.vb); %0.5
p.sa=1;
p.sb=1;
p.sc=1;
p.sv=1;


% Bounds for std hyperparameters (sig)
p.sigbd=[1e-3 1.4]; % 1.4

% Hyperprior location and scale for the prior spread hyperparameters (tau)
p.ea=p.glogit(0.7,p.sigbd); % Was 0.5. Corresponds to eta_alpha in the table
p.eb=p.glogit(0.7,p.sigbd); % 
p.ec=p.glogit(0.7,p.sigbd);
p.ev=p.glogit(0.7,p.sigbd);
p.xa=1; % Was 0.5 Corresponds to chi_alpha in the table
p.xb=1;
p.xc=1;
p.xv=1;

% Merge these into vectors for later convenience.
p.m=[p.ma; p.mb; p.mc; p.mv];
p.s=[p.sa; p.sb; p.sc; p.sv];
p.e=[p.ea; p.eb; p.ec; p.ev];
p.x=[p.xa; p.xb; p.xc; p.xv];




           
dovis=1;
p.vispri=0;

if dovis
%% Visualize the hyperpriors
fname='Helvetica';
tr=0.5;
rng(123);
close all;

lw=0.5;
%% b: temp bias
normal=@(x,mu,sd) (1/sqrt(2*pi*sd^2)).*exp(-0.5.*((x-mu)./sd).^2);
ncdf=@(x,mu,sd) 0.5.*(1+...
    erf((x-mu)./(sqrt(2.*sd.^2))));

Ns=1e1;

fs=16;
%figure('units','normalized','outerposition',[0 0 1 1])
fis=figure('units','inch','position',[0,0,18,8]);

tld = tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding', 'compact');
%subplot(2,4,3);
nexttile(2);
ppc=[15/255 82/255 186/255];
b=-5:0.05:5; b=b';
mu=zeros(Ns,1); sd=mu;
for j=1:Ns
    muj=p.mb+p.sb.*randn(1,1); mu(j)=muj;
    sdj=p.gexpit(p.eb+p.xb.*randn(1,1),p.sigbd); sd(j)=sdj;
    pj=normal(b,muj,sdj);
    plot(b,pj,'LineWidth',lw,'Color',[ppc tr]); hold on;
end
pj=normal(b,p.mb,p.gexpit(p.eb,p.sigbd));
plot(b,pj,'LineWidth',4*lw,'Color',[ppc 1]);

%tl='Hyperprior ensemble $p(\mu_b,\sigma_b)$';
%title(tl,'Interpreter','Latex','FontSize',fs);
tl='(b)';
title(tl,'Interpreter','tex','FontSize',fs);
xl='\textsf{Temperature bias, $b$ [K]}';
%xlabel(xl,'Interpreter','Latex','FontSize',fs);
%yl='PDF $p(b|\mu_\beta^{(i)},\tau_\beta^{(i)})$, [K]$^{-1}$';
%yl='Prior probability density';
%ylabel(yl,'Interpreter','Latex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01]);
set(groot, 'defaultAxesTickLabelInterpreter','tex');
%set(gca,'TickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
axis square;
%ylim([0 1]);
xlim([-4 4]);
%set(gca,'yscale','log');

%subplot(2,4,4);
nexttile(4+2);
cc=[0.4 0 0.8];
for j=1:Ns
    muj=mu(j); sdj=sd(j);
    cj=ncdf(b,muj,sdj);
    plot(b,cj,'LineWidth',lw,'Color',[cc tr]); hold on;
end
cj=ncdf(b,p.mb,p.gexpit(p.eb,p.sigbd));
plot(b,cj,'LineWidth',4*lw,'Color',cc);


%tl='Hyperprior ensemble $p(\mu_b,\sigma_b)$';
%title(tl,'Interpreter','Latex','FontSize',fs);
tl='(f)';
title(tl,'Interpreter','tex','FontSize',fs);
xl='Temperature bias [K]';
xts=-4:1:4;
set(gca,'XTick',xts);
xlabel(xl,'Interpreter','tex','FontSize',fs);
%yl='CDF $P(b|\mu_\beta^{(i)},\tau_\beta^{(i)})$';
%yl='Cumulative distribution function [-]';
%ylabel(yl,'Interpreter','Latex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01]);
set(groot, 'defaultAxesTickLabelInterpreter','tex');
%set(gca,'TickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
axis square;
%ylim([0 1]);
xlim([-4 4]);




%% a: Degree day factor


glogit=@(x,xbds) log((x-xbds(1))./diff(xbds))-log((xbds(2)-x)./diff(xbds));

glogitnormal=@(x,mu,sd,xbds) (1/sqrt(2*pi*sd^2)).*...
    (diff(xbds)./((x-xbds(1)).*(xbds(2)-x))).*...
    exp(-0.5.*((glogit(x,xbds)-mu)./sd).^2);

gltcdf=@(x,mu,sd,xbds) 0.5.*(1+...
    erf((glogit(x,xbds)-mu)./(sqrt(2.*sd.^2))));


%subplot(2,4,1);
nexttile(1);

a=p.ab(1):0.05:p.ab(2); a=a';
mu=zeros(Ns,1);
sd=zeros(Ns,1);
for j=1:Ns
    muj=p.ma+p.sa.*randn(1,1);
    mu(j)=muj;
    sdj=p.gexpit(p.ea+p.xa.*randn(1,1),p.sigbd);
    sd(j)=sdj;
    pj=glogitnormal(a,muj,sdj,p.ab);
    pj(a==p.ab(1)|a==p.ab(2))=0;
    %pj=gltcdf(a,muj,sdj,p.ab);
    plot(a,pj,'LineWidth',lw,'Color',[ppc tr]); hold on;
end
pj=glogitnormal(a,p.ma,p.gexpit(p.ea,p.sigbd),p.ab);
pj(a==p.ab(1)|a==p.ab(2))=0;
plot(a,pj,'LineWidth',4*lw,'Color',[ppc 1]); hold on;

%tl='Hyperprior ensemble $p(\mu_a,\sigma_a)$';
%title(tl,'Interpreter','Latex','FontSize',fs);
tl='(a)';
title(tl,'Interpreter','tex','FontSize',fs);
%xl='Degree day factor $a$ [mm K$^{-1}$ day$^{-1}$]';
%xlabel(xl,'Interpreter','Latex','FontSize',fs);
%yl='PDF $p(a|\mu_\alpha^{(i)},\tau_\alpha^{(i)})$, [mm K$^{-1}$ day$^{-1}$]$^{-1}$';
yl='Probability density';
ylabel(yl,'Interpreter','tex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01]);
set(groot, 'defaultAxesTickLabelInterpreter','tex');
%set(gca,'TickLabelInterpreter','latex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
axis square;
ylim([0 1]);
xlim(p.ab);


%subplot(2,4,2);
nexttile(4+1);
for j=1:Ns
    muj=mu(j);
    sdj=sd(j);
    cj=gltcdf(a,muj,sdj,p.ab);
    plot(a,cj,'LineWidth',lw,'Color',[cc tr]); hold on;
end
cj=gltcdf(a,p.ma,p.gexpit(p.ea,p.sigbd),p.ab);
plot(a,cj,'LineWidth',4*lw,'Color',[cc 1]); hold on;

%tl='Hyperprior ensemble $p(\mu_a,\scumuligma_a)$';
%title(tl,'Interpreter','Latex','FontSize',fs);
tl='(e)';
title(tl,'Interpreter','tex','FontSize',fs);
xl='Degree day factor [mm/(K day)]';
xlabel(xl,'Interpreter','tex','FontSize',fs);
yl='Cumulative probability';
ylabel(yl,'Interpreter','tex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01]);
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
%set(gca,'TickLabelInterpreter','latex');
axis square;
%ylim([0 1]);
xlim(p.ab);


%% c: Snowfall scaling

%subplot(2,4,5);
nexttile(3);

c=min(p.cb):0.025:max(p.cb); c=c';
mu=zeros(Ns,1); sd=mu;
for j=1:Ns
    muj=p.mc+p.sc.*randn(1,1); mu(j)=muj;
    if Ns==1
        muj=glogit(5,p.cb); mu(j)=muj;
    end
    sdj=p.gexpit(p.ec+p.xc.*randn(1,1),p.sigbd); sd(j)=sdj;
    pj=glogitnormal(c,muj,sdj,p.cb);
    pj(c==p.cb(1)|c==p.cb(2))=0;
    plot(c,pj,'LineWidth',lw,'Color',[ppc tr]); hold on;
end
pj=glogitnormal(c,p.mc,p.gexpit(p.ec,p.sigbd),p.cb);
pj(c==p.cb(1)|c==p.cb(2))=0;
plot(c,pj,'LineWidth',4*lw,'Color',[ppc 1]); hold on;

tl='(c)';
title(tl,'Interpreter','tex','FontSize',fs);
%xlabel(xl,'Interpreter','Latex','FontSize',fs);
%yl='PDF $p(c|\mu_\gamma^{(i)},\tau_\gamma^{(i)})$';
yl='Probability density';
%ylabel(yl,'Interpreter','Latex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01],...
    'YScale','linear');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
axis square;
ylim([0 3]);
xlim(p.cb);

%yyaxis right;

%subplot(2,4,6);
nexttile(4+3);
for j=1:Ns
    muj=mu(j);
    sdj=sd(j);
    cj=gltcdf(c,muj,sdj,p.cb);
    plot(c,cj,'LineWidth',lw,'Color',[cc tr]); hold on;
end
cj=gltcdf(c,p.mc,p.gexpit(p.ec,p.sigbd),p.cb);
plot(c,cj,'LineWidth',4*lw,'Color',[cc 1]); hold on;

%tl='Hyperprior ensemble $p(\mu_c,\sigma_c)$';
%title(tl,'Interpreter','Latex','FontSize',fs);
tl='(g)';
title(tl,'Interpreter','tex','FontSize',fs);
xl='Snowfall scaling [-]';
xlabel(xl,'Interpreter','tex','FontSize',fs);
%yl='CDF $P(c|\mu_\gamma^{(i)},\tau_\gamma^{(i)})$';
yl='Cumulative probability';
%ylabel(yl,'Interpreter','Latex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01],...
    'YScale','linear');
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
%set(gca,'TickLabelInterpreter','tex');
axis square;
ylim([0 1]);
xlim(p.cb);


%% v: CV 

%subplot(2,4,7);
nexttile(4);

v=0.01:0.01:1; v=v';
mu=zeros(Ns,1); sd=zeros(Ns,1);
for j=1:Ns
    muj=p.mv+p.sv.*randn(1,1); mu(j)=muj;
    sdj=p.gexpit(p.ev+p.xv.*randn(1,1),p.sigbd); sd(j)=sdj;
    pj=glogitnormal(v,muj,sdj,p.vb);
    pj(v==p.vb(1)|v==p.vb(2))=0;
    plot(v,pj,'LineWidth',lw,'Color',[ppc tr]); hold on;
end
pj=glogitnormal(v,p.mv,p.gexpit(p.ev,p.sigbd),p.vb);
pj(v==p.vb(1)|v==p.vb(2))=0;
plot(v,pj,'LineWidth',4*lw,'Color',[ppc 1]); hold on;

tl='(d)';
title(tl,'Interpreter','tex','FontSize',fs);
%xl='Ceofficient of variation {\it v} [-]';
%xlabel(xl,'Interpreter','Latex','FontSize',fs);
%yl='PDF $p(v|\mu_\nu^{(i)},\tau_\nu^{(i)})$';
%yl='Probability density';
%ylabel(yl,'Interpreter','Latex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01]);
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
axis square;
%ylim([0 2]);
xlim(p.vb);


%subplot(2,4,8);
nexttile(8);

for j=1:Ns
    muj=mu(j);
    sdj=sd(j);
    cj=gltcdf(v,muj,sdj,p.vb);
    plot(v,cj,'LineWidth',lw,'Color',[cc tr]); hold on;
end
cj=gltcdf(v,p.mv,p.gexpit(p.ev,p.sigbd),p.vb);
plot(v,cj,'LineWidth',4*lw,'Color',[cc 1]); hold on;



tl='(h)';
title(tl,'Interpreter','tex','FontSize',fs);
xl='Coefficient of variation [-]'; %xl='Coefficient of variation {\itv} [-]';
xlabel(xl,'Interpreter','tex','FontSize',fs);
%yl='CDF $P(v|\mu_\nu^{(i)},\tau_\nu^{(i)})$';
%yl='Cumulative probability';
%ylabel(yl,'Interpreter','Latex','FontSize',fs);
set(gca,'TickDir','out','LineWidth',2,'TickLength',[0.005, 0.01]);
set(groot, 'defaultAxesTickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex');
set(gca,'TickLabelInterpreter','tex','FontName',fname);
axis square;
%ylim([0 2]);
xlim(p.vb);

%%
%return

pause(2);
ax=gca;
ax.Toolbar.Visible = 'off';
toprint=sprintf('hyperpriors.pdf');
%exportgraphics(gcf,toprint,'Resolution',400);
exportgraphics(tld, toprint, 'ContentType', 'vector');

end
