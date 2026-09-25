clearvars;
close all;
tardir='results/';
% Check input data exists, otherwise fetch it from zenodo
if ~any(exist(tardir,'dir'))
    get_results;
end
set(0, 'DefaultFigureRenderer', 'painters');
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultAxesTitleFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');

fname='Helvetica';

methods={'CP'; 'NP'; 'PP'; 'MAP'};
ssites={'HS'; 'AK'; 'FF'; 'WF'};
gsites={'AB'; 'SG'; 'SB'; 'CD'};


load('~/work/ScientificColourMaps7/batlow/CategoricalPalettes/batlowS.mat');

tr=0.1;

fis=figure('units','inch','position',[0,0,18,8]);
pause(2);

c2p=1;
experiment='D0_F1_LO0';
%experiment='D1_F0_LO1';
compare={'cp','np','pp'};
figname=sprintf('experiment_%s',experiment);

fa=0.2;
cpc=[0.5 0.5 0.5];
cpm='s';
npc=[0.8 0 0]; 
npm='d';
ppc=[15/255 82/255 186/255];
ppm='c';
msz=40;
ec='none';

fssites={'Hornsund';'Abisko';'Filefjell';'Weissfluhjoch'};

pan='abcdefgh';


tld = tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding', 'compact');

for s=1:4
    site=ssites{s};
    fsite=fssites{s};
    tarf=sprintf('%s/Snow_%s_%s.mat',...
        tardir,site,experiment);
    load(tarf);



    obsfile=sprintf('input/snow/%s/%s_obs.mat',site,site);
    load(obsfile);
    disp('assigning obs');
    r.obs.t=obs.t; % Overwrite with non thinned obs
    r.obs.D=obs.D;



    if any(strcmp(site,'HS'))
        if any(strcmp(experiment(1:2),'D0'))
            ll=[0 800]; 
        else
            ll=[0 800]; 
        end
        ticks=0:100:ll(end); 
    elseif any(strcmp(site,'AK'))
        ticks=0:50:600; 
        if any(strcmp(experiment(1:2),'D0'))
            ll=[0 900]; 
        else
            ll=[0 800]; 
        end
        ticks=0:100:ll(end);
    elseif any(strcmp(site,'FF'))
        if any(strcmp(experiment(1:2),'D0'))
            ll=[0 1200]; 
        else
            ll=[0 1000];
        end
        ticks=0:100:ll(end);
    elseif any(strcmp(site,'WF'))
        if any(strcmp(experiment(1:2),'D0'))
            ll=[0 2000]; 
        else
            ll=[0 1600]; 
        end
        ticks=0:200:ll(end);
    end

    cpce=[0 0 0];
    npce=[0 0 0];
    ppce=[0 0 0];

    t=datenum(sprintf('01-Sep-%d',r.p.wys(1)-1)):...
        (datenum(sprintf('01-Sep-%d',r.p.wys(end)))-1);
    t=t';
    Nt=numel(t);
    twy=datevec(t); 
    next=twy(:,2)>9; twy=twy(:,1); twy(next,1)=twy(next,1)+1;

    % Assume all obs dates exist
    % This approach also handles cases with multiple obs per day
    No=numel(r.obs.t);
    theseo=zeros(No,1); % indices not logical array allows for copies
    tinds=(1:Nt)';
    for j=1:No
        herej=t==r.obs.t(j);
        theseo(j)=tinds(herej);
    end
    wyo=datevec(r.obs.t);
    next=wyo(:,2)>=9; wyo=wyo(:,1); wyo(next,1)=wyo(next,1)+1;
    wyu=unique(wyo);
    if any(strcmp(experiment,'D0_F1_LO0'))
        theseoa=wyo>=2000&wyo<=2023;
    else
        theseoa=wyo<r.p.lowys(1);
    end

    
    fs=18;
    if s==1
        nexttile;
    else
        nexttile(s);
    end
    x=r.obs.D;
    if any(strcmp(compare,'cp'))
        y1=r.cp.Dm(theseo); y1s=r.cp.Ds(theseo);
        y2=r.np.Dm(theseo); y2s=r.np.Ds(theseo);
    end
    y3=r.pp.Dm(theseo); y3s=r.pp.Ds(theseo);
    Ns=1e3;
    y0=r.pp.Dprim(theseo); y0s=r.pp.Dpris(theseo);
    y0d=max(y0+1.*y0s.*randn(size(y0,1),Ns),0);
    y1d=max(y1+1.*y1s.*randn(size(y1,1),Ns),0);
    y2d=max(y2+1.*y2s.*randn(size(y2,1),Ns),0);
    y3d=max(y3+1.*y3s.*randn(size(y3,1),Ns),0);
    xx=max([max(x) max(y1) max(y2)]);
    nn=min([min(x) min(y1) min(y2)]);
    xlim(ll); ylim(ll);
    fa=1;
    hold on;
    plot(ll,ll,'-k');


    % Binning procedure
    xp=x(theseoa); 
    xpv=x(~theseoa);
    xx=max(x)+eps;
    if s<0
        dx=10; % bin size
    elseif s<4
        dx=25;
    else
        dx=50;
    end
    xx=ceil(xx/dx)*dx;
    nb=xx/dx;

    xbv=nan(nb,1);
    ymc=nan(nb,4);
    ysc=nan(nb,4);
    ymv=nan(nb,4);
    ysv=nan(nb,4);

    % Add a zero bin
    nb=nb+1;
    for b=1:nb
        if b==1
            bmin=0; 
            bmax=1e-6;
        else
            bmin=bmax;
            bmax=(b-1)*dx;
        end
        theseb=xp>=bmin&xp<bmax;
        thesebv=xpv>=bmin&xpv<bmax;
        if b==1
            xbm=0;
        else
            xbm=bmin+dx/2;
        end
        xbv(b)=xbm;

        % Order: pri,CP, NP, PP
        for approach=0:3
            if approach==0
                yp=y0d(theseoa,:);
                ypv=y0d(~theseoa,:);
            elseif approach==1
                yp=y1d(theseoa,:);
                ypv=y1d(~theseoa,:);
            elseif approach==2
                yp=y2d(theseoa,:);
                ypv=y2d(~theseoa,:);
            else
                yp=y3d(theseoa,:);
                ypv=y3d(~theseoa,:);
            end
            yp=yp(theseb,:);
            ypv=ypv(thesebv,:);
            if isempty(yp)||sum(theseb)<5
                yp=nan;
            end
            if isempty(ypv)||sum(thesebv)<5
                ypv=nan;
            end
            ymc(b,approach+1)=mean(yp(:));
            ysc(b,approach+1)=std(yp(:));
            ymv(b,approach+1)=mean(ypv(:));
            ysv(b,approach+1)=std(ypv(:));
        end

    end

    prc=[0 0 0];
    cpc=[0.5 0.5 0.5];
    npc=[0.800, 0.475, 0.655]; 
    ppc=[0 0 0.8]; 

    cs=[prc;...
        cpc;...
        npc;...
        ppc];
    alp=0.2;
    Napp=4;


    aorder=[1 2 3 4]; % Order: From wide to narrow.
    for approach=1:Napp
        ao=aorder(approach);
        ca=cs(ao,:);
        ym=ymc(:,ao);
        ys=ysc(:,ao);
        cp=[ym-ys ym ym+ys];
        these=~isnan(ym);
        xb=xbv(these);
        cp=cp(these,:);
        % Skip gaps
        here=diff(xb)>dx;
        bpoint=xb(here);
        bpoint=[bpoint; max(xb)];
        bpoint=unique(bpoint);
        bjold=-Inf;
        for j=1:numel(bpoint)
            bj=bpoint(j);
            thesej=xb>bjold&xb<=bj;
            if ao>1
                alpj=alp;
                caf=ca;
                lst='-';
            else
                alpj=alp*0.5;
                caf=[1 1 1];
                lst=':';
            end
            percp(ao)=percplot(xb(thesej)',cp(thesej,:)',caf,alpj,lst,ca); 
            bjold=bj;
            if ao==1
                tmpp(j)=percp(ao);
            end
        end

    end

    % Validation order should be different
    aorder=[1 2 3 4];
    for approach=1:Napp
        ao=aorder(approach);
        if ao>1
            lst='-';
            lw=2;
            alpis=2*alp;
        else
            lst='--';
            lw=1;
            alpis=1;
        end
        if ao>0
            plot(xbv,ymc(:,ao),'Color',[cs(ao,:) alpis],'LineWidth',lw,...
                'LineStyle',lst);
        end
    end

     

    xlim(ll); ylim(ll);
    axis square;
    box on; 
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','tex');
    set(gca,'TickLabelInterpreter','tex','FontName',fname);
    set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
    set(gca,'XTick',ticks,'YTick',ticks);
    if any(strcmp(experiment,'D0_F1_LO0'))
        calstr='calibration';
        expstr='Exp. F';
    else
        calstr='calibration';
        expstr='Exp. S,';
    end
    tis=sprintf('(%s) %s %s',pan(s),fsite,calstr);
    title(tis,'Interpreter','tex','FontSize',0.8.*fs);

    if s==1
        leg=legend(percp,{'Prior';'Complete pooling';'No pooling';'Partial pooling'},...
            'Interpreter','tex','FontSize',12,...
            'Location','NorthEast');
    end
    
    nexttile(4+s);
    hold on;
    plot(ll,ll,'-k');

     % Calibration order should be different
   
     aorder=[1 3 4 2]; % Order: From wide to narrow.
     for approach=1:Napp
        ao=aorder(approach);
        ca=cs(ao,:);
        ym=ymv(:,ao);
        ys=ysv(:,ao);
        cp=[ym-ys ym ym+ys];
        these=~isnan(ym);
        xb=xbv(these);
        cp=cp(these,:);
        % Skip gaps
        here=diff(xb)>dx;
        bpoint=xb(here);
        bpoint=[bpoint; max(xb)];
        bpoint=unique(bpoint);
        bjold=-Inf;
        for j=1:numel(bpoint)
            bj=bpoint(j);
            thesej=xb>bjold&xb<=bj;
            if ao>1
                alpj=alp;
                caf=ca;
                lst='-';
                lw=2;
            else
                alpj=alp*0.5;
                caf=[1 1 1];
                lst=':';
                lw=1;
            end

            percp(ao)=percplot(xb(thesej)',cp(thesej,:)',caf,alpj,lst,ca); 
            bjold=bj;
            if ao==1
                tmp(j)=percp(ao);
            end
        end

     end


    for approach=1:4
        ao=aorder(approach);
        if ao>1
            lst='-';
            lw=2;
            alpis=2*alp;
        else
            lst='--';
            lw=1;
            alpis=1;
        end
        if ao>0
            plot(xbv,ymv(:,ao),'Color',[cs(ao,:) alpis],'LineWidth',lw,...
                'LineStyle',lst);
        end
    end

      
    

    xlim(ll); ylim(ll);
    axis square;
    box on; 
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','tex');
    set(gca,'TickLabelInterpreter','tex','FontName',fname);
    set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
    set(gca,'XTick',ticks,'YTick',ticks);
    if any(strcmp(experiment,'D0_F1_LO0'))
        calstr='validation';
    else
        calstr='validation';
    end
    tis=sprintf('(%s) %s %s',pan(4+s),fsite,calstr);
    title(tis,'Interpreter','tex','FontSize',0.8.*fs);



end

xlabel(tld, 'Observed snow water equivalent [mm w.e.]', 'Interpreter','tex','FontSize',1.2*fs);
ylabel(tld, 'Estimated snow water equivalent [mm w.e.]', 'Interpreter','tex','FontSize',1.2*fs);
if any(strcmp(experiment,'D0_F1_LO0'))
    titleis='Experiment F: FSCA assimilation';
else
    titleis='Experiment S: SWE assimilation';
end


ax=gca;
ax.Toolbar.Visible = 'off';
pause(2);
toprint=sprintf('results/%s.pdf',figname);
exportgraphics(tld, toprint, 'ContentType', 'vector');




%% Glacier scatter


tr=0.2;
close all;
fis=figure('units','inch','position',[0,0,18,8]);
fgsites={'Austre Brøggerbreen';'Storglaciären';'Storbrean';'Claridenfirn'};



tld = tiledlayout(2,4, 'TileSpacing', 'compact', 'Padding', 'compact');



pause(2);
experiment='LO1';
for s=1:4
    site=gsites{s};
    fsite=fgsites{s};
    tarf=sprintf('%s/Glacier_%s_%s.mat',...
        tardir,site,experiment);
    load(tarf);

    t=datenum(sprintf('01-Sep-%d',r.p.wys(1)-1)):...
        (datenum(sprintf('01-Sep-%d',r.p.wys(end)))-1);
    t=t';
    Nt=numel(t);
    twy=datevec(t); 
    next=twy(:,2)>9; twy=twy(:,1); twy(next,1)=twy(next,1)+1;

    % Assume all obs dates exist
    % This approach also handles cases with multiple obs per day
    No=numel(r.obs.t);
    theseo=zeros(No,1); % indices not logical array allows for copies
    tinds=(1:Nt)';
    for j=1:No
        herej=t==r.obs.t(j);
        theseo(j)=tinds(herej);
    end
    wyo=datevec(r.obs.t);
    next=wyo(:,2)>9; wyo=wyo(:,1); wyo(next,1)=wyo(next,1)+1;
    theseoa=wyo<r.p.lowys(1);

    prc=[0.3 0.3 0.3];
    cpce=cpc;
    npce=npc;
    ppce=ppc;
    prce=prc;

    if any(strcmp(site,'AB'))
        
        ll=[-3000 1500];
        ticks=ll(1):500:ll(2);
    elseif any(strcmp(site,'SG'))
        ll=[-5000 3000]; 
        ticks=ll(1):1000:ll(2);
    elseif any(strcmp(site,'SB'))
        ll=[-7000 3000]; 
        ticks=ll(1):1000:ll(2);
    elseif any(strcmp(site,'CD'))
        ll=[-5000 3000];
        ticks=ll(1):1000:ll(2);
    end

    
    fs=18;
    nexttile(s);
    
    x=r.obs.D;
    y0=r.pp.Dprim(theseo); y0s=r.pp.Dpris(theseo);
    y1=r.cp.Dm(theseo); y1s=r.cp.Ds(theseo);
    y3=r.np.Dm(theseo); y3s=r.np.Ds(theseo);
    y2=r.pp.Dm(theseo); y2s=r.pp.Ds(theseo);
    xx=max([max(x) max(y1) max(y2)]);
    nn=min([min(x) min(y1) min(y2)]);

    trp=0.4;
    xlim(ll); ylim(ll);
    err1=y1(theseoa)-x(theseoa);
    rmse1=sqrt(mean(err1.^2));
    r1=min(min(corrcoef(y1(theseoa),x(theseoa))));
    crps1=mean(CRPSg(y1(theseoa),y1s(theseoa),x(theseoa)));
    err2=y2(theseoa)-x(theseoa);
    rmse2=sqrt(mean(err2.^2));
    crps2=mean(CRPSg(y2(theseoa),y2s(theseoa),x(theseoa)));
    r2=min(min(corrcoef(y2(theseoa),x(theseoa))));


    Ns=5e0;
    y1d=y1+y1s.*randn(size(y1,1),Ns);
    y2d=y2+y2s.*randn(size(y2,1),Ns);

    plot([0 0],ll,'-k','LineWidth',0.25); hold on;
    plot(ll,[0 0],'-k','LineWidth',0.25);
    plot(ll,ll,'-k'); 

    plot([x(theseoa)'; x(theseoa)'],...
        [(y0(theseoa)+y0s(theseoa))'; (y0(theseoa)-y0s(theseoa))'],...
        'LineStyle','-','LineWidth',0.25,'Color',[prc 0.1]); hold on;
    plot([x(theseoa)'; x(theseoa)'],...
        [(y1(theseoa)+y1s(theseoa))'; (y1(theseoa)-y1s(theseoa))'],...
        'LineWidth',0.5,'Color',[cpc trp]); hold on;
    plot([x(theseoa)'; x(theseoa)'],...
        [(y2(theseoa)-y2s(theseoa))' ; (y2(theseoa)+y2s(theseoa))'],...
        'LineWidth',1,'Color',[ppc trp]);
    plot([x(theseoa)'; x(theseoa)'],...
        [(y3(theseoa)-y3s(theseoa))' ; (y3(theseoa)+y3s(theseoa))'],...
        'LineWidth',0.75,'Color',[npc trp]);

    sc(1)=scatter(x(theseoa),y0(theseoa,:),msz,'hexagram',...
        'MarkerFaceColor',[1 1 1],'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',prc,'LineWidth',0.5); hold on;
    sc(2)=scatter(x(theseoa),y1(theseoa,:),1.5.*msz,cpm,...
        'MarkerFaceColor',cpc,'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',cpce,'LineWidth',0.5); hold on;
    sc(3)=scatter(x(theseoa),y3(theseoa,:),1.5.*msz,npm,...
        'MarkerFaceColor',npc,'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',npce,'LineWidth',0.5); hold on;
    sc(4)=scatter(x(theseoa),y2(theseoa,:),msz,ppm,...
        'MarkerFaceColor',ppc,'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',ppce,'LineWidth',0.5); hold on;
    xlim(ll); ylim(ll);
    axis square;
    set(gca,'XTick',ticks,'YTick',ticks);
    box on; 
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','tex');
    set(gca,'TickLabelInterpreter','tex','FontName',fname);
    set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
    tis=sprintf('(%s) %s calibration',pan(s),fsite);
    title(tis,'Interpreter','tex','FontSize',0.8.*fs);
    if s==4
        leg=legend(sc,{'Prior',;'Complete pooling';'No pooling';'Partial pooling'},...
            'Interpreter','tex','FontSize',12,...
            'Location','SouthEast');
    end
    
    

    nexttile(4+s);
    plot([0 0],ll,'-k','LineWidth',0.25); hold on;
    plot(ll,[0 0],'-k','LineWidth',0.25);

    plot(ll,ll,'-k'); hold on;

    plot([x(~theseoa)'; x(~theseoa)'],...
        [(y0(~theseoa)+y0s(~theseoa))'; (y0(~theseoa)-y0s(~theseoa))'],...
        'LineWidth',0.1,'Color',[prc 0.1]); hold on;


    plot([x(~theseoa)'; x(~theseoa)'],...
        [(y2(~theseoa)-y2s(~theseoa))' ; (y2(~theseoa)+y2s(~theseoa))'],...
        'LineWidth',1,'Color',[ppc trp]);
    plot([x(~theseoa)'; x(~theseoa)'],...
        [(y3(~theseoa)-y3s(~theseoa))' ; (y3(~theseoa)+y3s(~theseoa))'],...
        'LineWidth',0.75,'Color',[npc trp]);

    plot([x(~theseoa)'; x(~theseoa)'],...
        [(y1(~theseoa)+y1s(~theseoa))'; (y1(~theseoa)-y1s(~theseoa))'],...
        'LineWidth',0.5,'Color',[cpc trp]); hold on;
    

    scatter(x(~theseoa),y0(~theseoa,:),msz,'hexagram',...
        'MarkerFaceColor',[1 1 1],'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',prc,'LineWidth',0.5); hold on;

    scatter(x(~theseoa),y3(~theseoa,:),1.5.*msz,npm,...
        'MarkerFaceColor',npc,'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',npce,'LineWidth',0.5); hold on;
    scatter(x(~theseoa),y1(~theseoa,:),1.5.*msz,cpm,...
        'MarkerFaceColor',cpc,'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',cpce,'LineWidth',0.5); hold on;
    scatter(x(~theseoa),y2(~theseoa,:),msz,'o',...
        'MarkerFaceColor',ppc,'MarkerFaceAlpha',tr,...
        'MarkerEdgeColor',ppce,'LineWidth',0.5); hold on;
    xlim(ll); ylim(ll);
    err1=y1(~theseoa)-x(~theseoa);
    rmse1=sqrt(mean(err1.^2));
    crps1=mean(CRPSg(y1(~theseoa),y1s(~theseoa),x(~theseoa)));
    r1=min(min(corrcoef(y1(~theseoa),x(~theseoa))));
    err2=y2(~theseoa)-x(~theseoa);
    rmse2=sqrt(mean(err2.^2));
    crps2=mean(CRPSg(y2(~theseoa),y2s(~theseoa),x(~theseoa)));
    r2=min(min(corrcoef(y2(~theseoa),x(~theseoa))));
    axis square;
    set(gca,'XTick',ticks,'YTick',ticks);
    box on; 
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.005, 0.005]);
    set(groot, 'defaultAxesTickLabelInterpreter','tex');
    set(gca,'TickLabelInterpreter','tex','FontName',fname);
    set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
  
    tis=sprintf('(%s) %s validation',pan(4+s),fsite);
    title(tis,'Interpreter','tex','FontSize',0.8.*fs);
  

end


xlabel(tld, 'Observed glacier-wide mass balance [mm w.e.]', 'Interpreter','tex','FontSize',1.2*fs);
ylabel(tld, 'Estimated glacier-wide mass balance [mm w.e.]', 'Interpreter','tex','FontSize',1.2*fs);
ax=gca;
ax.Toolbar.Visible = 'off';
pause(2);

figname=sprintf('glacier_scatter');

toprint=sprintf('results/%s.pdf',figname);
exportgraphics(tld, toprint, 'ContentType', 'vector');