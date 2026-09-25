clearvars;
tardir='results/';

% Check input data exists, otherwise fetch it from zenodo
if ~any(exist(tardir,'dir'))
    get_results;
end
close all;
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultAxesTitleFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
fis=figure('units','inch','position',[0,0,18,18]);
pause(2);
pans='abcdefgh';
pars={'Degree day factor [mm K^{-1} day^{-1}]';...
    'Air temperature bias [K]';...
    'Snowfall scaling [-]'};
parbs=[1 7;
       -2 4;
       0 5];
p=3;
fname='Helvetica';

doglac=1;
if doglac
    tld = tiledlayout(4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
else
    tld = tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
end
for j=1:4
    sites={'HS';'AK';'FF';'WF'};
    fsites={'Hornsund';'Abisko';'Filefjell';'Weissfluhjoch'};
   
    yl=parbs(p,:);
    if p==3
        yl(2)=3.5;
    end
    
    nexttile;
    
    
    fa=1;
    prc=[0 0 0];
    cpc=[0.5 0.5 0.5];
    npc=[0.800 0.475 0.655]; 
    ppc=[0 0 0.8]; 
    ec='none';

    nsig=1;
    alp=0.2;

    plot([2000 2000],yl,'Color',[0.5 0.5 0.5 0.4],'LineWidth',0.5,'LineStyle','-'); hold on;
    load(sprintf([tardir 'Snow_%s_D1_F0_LO0.mat'],sites{j}));
    hold on;
    t=r.p.wys;
    these=t>=2000;
    refm=r.pp.phim(:,p);
    refm=refm(these);
    refs=r.pp.phis(:,p);
    refs=refs(these);
    tref=t(these);
    ref=[refm-nsig.*refs refm refm+nsig.*refs];
    refc=[0.95 0.8 0.2];

    load(sprintf([tardir 'Snow_%s_D1_F0_LO1.mat'],sites{j}));


    % Assuming bounded, modify for b
    Nh=1e3; % Number of hyperparameter set combinations to sample
    Ne=1e3; % Number of parameters to sample for each hyperparameter set
    mu=r.p.m(p)+r.p.s(p).*randn(Nh,1);
    sig=r.p.e(p)+r.p.x(p).*randn(Nh,1);
    prp=mu+sig.*randn(Nh,Ne);
    if p~=2
        if p==1
            bds=r.p.ab;
        elseif p==3
            bds=r.p.cb;
        end
        prp=r.p.gexpit(prp,bds);
    end
    prpm=mean(prp(:));
    prps=std(prp(:));
    clear prp;
    prpm=repmat(prpm,1,size(r.np.phim,1))';
    prps=repmat(prps,1,size(r.np.phim,1))';
    pr=[prpm-nsig.*prps prpm prpm+nsig.*prps];
    percp(1)=percplot(t,pr,[1 1 1],alp,':',prc); 

    cpm=repmat(r.cp.phim(p),1,size(r.np.phim,1),1)';
    cps=repmat(r.cp.phis(p),1,size(r.np.phim,1),1)';
    cp=[cpm-nsig.*cps cpm cpm+nsig.*cps];
    
    percp(2)=percplot(t,cp,cpc,alp,'-',cpc); 

    npm=r.np.phim(:,p);
    nps=r.np.phis(:,p);
    np=[npm-nsig.*nps npm npm+nsig.*nps]; hold on;

    percp(3)=percplot(t,np,npc,alp,'-',npc);  

    ppm=r.pp.phim(:,p);
    pps=r.pp.phis(:,p);
    pp=[ppm-nsig.*pps ppm ppm+nsig.*pps];
    percp(4)=percplot(t,pp,ppc,alp,'-',ppc); 


    percp(5)=percplot(tref,ref,refc,alp,'-',refc); hold on;

    plot(t,prpm,'--','LineWidth',1,'Color',[prc 1]); % 1
    pt(2)=plot(t,npm,'Color',[npc 2*alp],'LineWidth',2);
    pt(3)=plot(t,ppm,'Color',[ppc 2*alp],'LineWidth',2);
    pt(1)=plot(t,cpm,'Color',[cpc 2*alp],'LineWidth',2);

    pt(4)=plot(tref,refm,'Color',refc,'LineWidth',2);

    % Check coverage:
    here=t>=min(tref)&t<=max(tref);
    ppenv=pp(here,1:2:3);
    inside=refm>=ppenv(:,1)&refm<=ppenv(:,2);
    scoverage(j)=sum(inside)/numel(inside);

    npenv=np(here,1:2:3);
    inside=refm>=npenv(:,1)&refm<=npenv(:,2);
    npscoverage(j)=sum(inside)/numel(inside);

    cpenv=cp(here,1:2:3);
    inside=refm>=cpenv(:,1)&refm<=cpenv(:,2);
    cpscoverage(j)=sum(inside)/numel(inside);

    if p==1
        yts=1:7;
    elseif p==2
        yts=-6:1:6;
    elseif p==3
        yts=0:1:5;
    end
    set(gca,'YTick',yts);


    fs=18;
    box on; 
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.0025, 0.0025]);
    set(groot, 'defaultAxesTickLabelInterpreter','tex');
    set(gca,'TickLabelInterpreter','tex','FontName',fname);
    set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
    tis=sprintf('(%s) %s',pans((j-1)*2+1),fsites{j});
    title(tis,'Interpreter','tex','FontSize',1.*fs);
    xl=[min(t)-1 max(t)+1];
    xlim(xl);
    if j==0
        leg=legend(percp,{'Prior','Complete pooling','No pooling',...
            'Partial pooling','Reference'},...
            'Interpreter','tex','FontSize',0.7*fs,...
            'Location','NorthWest','Orientation','horizontal');
    end
    if p==2
       ylim([-6 3.5]); 
    else
        ylim(yl);
    end

    sites={'AB';'SG';'SB';'CD'};
    fsites={'Austre Brøggerbreen';'Storglaciären';'Storbrean';'Claridenfirn'};
    yl=parbs(p,:);

    nexttile;
    load(sprintf('~/Downloads/resnew/Glacier_%s_LO0.mat',sites{j}));
    t=r.p.wys;
    refm=r.pp.phim(:,p);
    these=t>=2000;
    refm=refm(these);
    refs=r.pp.phis(:,p);
    refs=refs(these);
    tref=t(these);
    ref=[refm-nsig.*refs refm refm+nsig.*refs];
    load(sprintf('~/Downloads/resnew/Glacier_%s_LO1.mat',sites{j}));
    
    fa=1;
    ec='none';

    nsig=1;
    alp=0.2;

    plot([2000 2000],yl,'Color',[0.5 0.5 0.5 0.4],'LineWidth',0.5,'LineStyle','-'); hold on;


    % Assuming bounded, modify for b
    Nh=1e3; % Number of hyperparameter set combinations to sample
    Ne=1e3; % Number of parameters to sample for each hyperparameter set
    mu=r.p.m(p)+r.p.s(p).*randn(Nh,1);
    sig=r.p.e(p)+r.p.x(p).*randn(Nh,1);
    prp=mu+sig.*randn(Nh,Ne);
    if p~=2
        if p==1
            bds=r.p.ab;
        elseif p==3
            bds=r.p.cb;
        end
        prp=r.p.gexpit(prp,bds);
    end
    prpm=mean(prp(:));
    prps=std(prp(:));
    clear prp;
    prpm=repmat(prpm,1,size(r.np.phim,1))';
    prps=repmat(prps,1,size(r.np.phim,1))';
    pr=[prpm-nsig.*prps prpm prpm+nsig.*prps];
    percp(1)=percplot(t,pr,[1 1 1],alp,':',prc); 

    cpm=repmat(r.cp.phim(p),1,size(r.np.phim,1),1)';
    cps=repmat(r.cp.phis(p),1,size(r.np.phim,1),1)';
    cp=[cpm-nsig.*cps cpm cpm+nsig.*cps];
    
    percp(2)=percplot(t,cp,cpc,alp,'-',cpc); 

    npm=r.np.phim(:,p);
    nps=r.np.phis(:,p);
    np=[npm-nsig.*nps npm npm+nsig.*nps]; hold on;

    percp(3)=percplot(t,np,npc,alp,'-',npc);  

    ppm=r.pp.phim(:,p);
    pps=r.pp.phis(:,p);
    pp=[ppm-nsig.*pps ppm ppm+nsig.*pps];
    percp(4)=percplot(t,pp,ppc,alp,'-',ppc); 


    

    percplot(tref,ref,refc,alp,'-',refc); hold on;

    plot(t,prpm,'--','LineWidth',1,'Color',[prc 1]); % 1
    pt(2)=plot(t,npm,'Color',[npc 2*alp],'LineWidth',2);
    pt(3)=plot(t,ppm,'Color',[ppc 2*alp],'LineWidth',2);
    pt(1)=plot(t,cpm,'Color',[cpc 2*alp],'LineWidth',2);

    
    pt(4)=plot(tref,refm,'Color',refc,'LineWidth',2);

    % Check coverage:
    here=t>=min(tref)&t<=max(tref);
    ppenv=pp(here,1:2:3);
    inside=refm>=ppenv(:,1)&refm<=ppenv(:,2);
    gcoverage(j)=sum(inside)/numel(inside);

    npenv=np(here,1:2:3);
    inside=refm>=npenv(:,1)&refm<=npenv(:,2);
    npgcoverage(j)=sum(inside)/numel(inside);

    cpenv=cp(here,1:2:3);
    inside=refm>=cpenv(:,1)&refm<=cpenv(:,2);
    cpgcoverage(j)=sum(inside)/numel(inside);
    
    
    if j==1
        leg=legend(percp,{'Prior','Complete pooling','No pooling',...
            'Partial pooling','Reference'},...
            'Interpreter','tex','FontSize',0.7*fs,...
            'Location','NorthWest','Orientation','horizontal');
    end


    box on; 

    set(gca,'YTick',yts);
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.0025, 0.0025]);
    set(groot, 'defaultAxesTickLabelInterpreter','tex');
    set(gca,'TickLabelInterpreter','tex','FontName',fname);
    set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
    tis=sprintf('(%s) %s',pans((j-1)*2+2),fsites{j});
    title(tis,'Interpreter','tex','FontSize',1.*fs);

    xl=[min(t)-1 max(t)+1];
    xlim(xl);
    
    if p==2
        yl(2)=3.5;
    end
    ylim(yl);
end



ax=gca;
ax.Toolbar.Visible = 'off';

ylabel(tld,sprintf('%s',pars{p}),'Interpreter','tex','FontSize',1.5*fs,...
    'FontName',fname);

pause(2);
figname='pardyn';
toprint=sprintf('results/%s_par%d.pdf',figname,p);
exportgraphics(tld, toprint, 'ContentType', 'vector');