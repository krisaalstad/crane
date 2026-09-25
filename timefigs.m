clearvars; 
dolong=0;
if dolong
    dt=7;
else
    dt=3;
end
set(groot, 'defaultAxesFontWeight', 'normal');
set(groot, 'defaultAxesTitleFontWeight', 'normal');
set(groot, 'defaultTextFontWeight', 'normal');
tardir='results/';
% Check input data exists, otherwise fetch it from zenodo
if ~any(exist(tardir,'dir'))
    get_results;
end


fname='Helvetica';

pans='abcdefgh';
parbs=[0 10;
       -3 3;
       0 5];
p=1;

methods={'static inference via complete pooling';'annual inference via no pooling';'hierarchical inference via partial pooling'};


pric=[0 0 0];
posc=[0 0.3 0.8];
obsc=[0.95 0.8 0.2];
mc=[0.5 0.5 0.5];
alp=0.2;
ec='none';
alphabet = 'abc';
nsig=1;
fs=18;
for expt=1:3
    if expt==1
        expis=[tardir 'Snow_%s_D1_F0_LO1.mat'];
        expc='S';
        nsig=1;
    elseif expt==2
        expis=[tardir 'Snow_%s_D0_F1_LO0.mat'];
        expc='F';
        nsig=1;
    elseif expt==3
        expis=[tardir 'Glacier_%s_LO1.mat'];
        expc='G';
        nsig=1;
    end

    for s=1:4

        if expt<3
            sites={'HS';'AK';'FF';'WF'};
            fsites={'Hornsund';'Abisko';'Filefjell';'Weissfluhjoch'};
        else
            sites={'AB';'SG';'SB';'CD'};
            fsites={'Austre Brøggerbreen';'Storglaciären';'Storbrean';'Claridenfirn'};
        end

if expt~=2
    nplots=1;
else
    nplots=2;
end

load(sprintf(expis,sites{s}));
t=datenum(sprintf('01-Sep-%d',r.p.wys(1)-1)):...
    (datenum(sprintf('01-Sep-%d',r.p.wys(end)))-1);
t=t';
Nt=numel(t);
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
site=sites{s};

for pltis=1:nplots

close all;
fis=figure('units','inch','position',[0,0,18,18]);
pause(2);



tld = tiledlayout(3,1, 'TileSpacing', 'compact', 'Padding', 'compact');
for j=1:3
    nexttile;

    plotF=expt==2&&pltis==2;

    if plotF
        to=r.obs.ts;
        Do=r.obs.Fs;
        lw=0.5;
    else
        to=r.obs.t;
        Do=r.obs.D;
        lw=0.5;
    end
    if j==1
        Dprm=r.cp.Dprim;
        Dprs=r.cp.Dpris;
        Dpom=r.cp.Dm;
        Dpos=r.cp.Ds;
        if plotF
            Dprm=r.cp.Fprim;
            Dprs=r.cp.Fpris;
            Dpom=r.cp.Fm;
            Dpos=r.cp.Fs;
        else
            Dprm=r.cp.Dprim;
            Dprs=r.cp.Dpris;
            Dpom=r.cp.Dm;
            Dpos=r.cp.Ds;
        end
    elseif j==2
        if plotF
            Dprm=r.np.Fprim;
            Dprs=r.np.Fpris;
            Dpom=r.np.Fm;
            Dpos=r.np.Fs;
        else
            Dprm=r.np.Dprim;
            Dprs=r.np.Dpris;
            Dpom=r.np.Dm;
            Dpos=r.np.Ds;
        end
    elseif j==3
        if plotF
            Dprm=r.pp.Fprim;
            Dprs=r.pp.Fpris;
            Dpom=r.pp.Fm;
            Dpos=r.pp.Fs;
        else
            Dprm=r.pp.Dprim;
            Dprs=r.pp.Dpris;
            Dpom=r.pp.Dm;
            Dpos=r.pp.Ds;
        end
    end

    if dt>1&&expt<3
        Dpom=movmean(Dpom,dt,1);
        Dprm=movmean(Dprm,dt,1);
        Dpos=movmean(Dpos,dt,1);
        Dprs=movmean(Dprs,dt,1);
    end

  

   
    Dpri=[Dprm-nsig.*Dprs Dprm Dprm+nsig.*Dprs];
    Dpst=[Dpom-nsig.*Dpos Dpom Dpom+nsig.*Dpos];
    if expt<3
        Dpri=max(Dpri,0);
        Dpst=max(Dpst,0);
    end
    
    if plotF
        Dpri=min(Dpri,1);
        Dpst=min(Dpst,1);
    end

     if any(strcmp(sites{s},'HS'))||any(strcmp(sites{s},'AK'))
        yls=[0 800]; 
    elseif any(strcmp(sites{s},'FF'))
        yls=[0 1000];
    elseif any(strcmp(sites{s},'WF'))
        yls=[0 1600];
    elseif any(strcmp(sites{s},'AB'))
         yls=[-2100 1100];
    elseif any(strcmp(sites{s},'SG'))
         yls=[-3100 3100];
    elseif any(strcmp(sites{s},'SB'))
         yls=[-5000 3500];
    elseif any(strcmp(sites{s},'CD'))
         yls=[-3500 3500];
    else
        yx=1.1.*max(abs(r.obs.D));
        yls=[-yx yx];
    end
    if plotF
        yls=[0 1];
    end
    ylim(yls);

    hold on;
    tvc=datenum('01-Sep-1999');
    plot([min(t) max(t)],[0 0],'Color',[0.5 0.5 0.5 0.4],'LineWidth',0.5);
    plot([tvc tvc],yls,'Color',[0.5 0.5 0.5 0.4],'LineWidth',0.5,'LineStyle','-'); hold on;
    


    if ~plotF
        price=pric;
        posce='none';
    else
        price=pric;
        posce='none';
    end
    if expt<3
        tmp=percplot(t(1:dt:end),Dpri(1:dt:end,:),pric,0.1*alp,':',price); hold on;
        percplot(t(1:dt:end),Dpst(1:dt:end,:),posc,alp,'-',posce);
        tmp.FaceAlpha=0.05;
        
    else
        Ny=numel(r.p.wys);
        for wyr=1:Ny
            tstart=datenum(sprintf('%s-%d',r.p.dstart,r.p.wys(wyr)-1));
            tend=datenum(sprintf('%s-%d',r.p.dend,r.p.wys(wyr)))-1;
            here=t>=tstart&t<=tend;
            hold on;
            th=t(here); Dprih=Dpri(here,:); Dpoh=Dpst(here,:);
            Dprmh=Dprm(here); Dpomh=Dpom(here);
            tmp(wyr)=percplot(th(1:dt:end),Dprih(1:dt:end,:),pric,0.1*alp,':',price); hold on;
            percplot(th(1:dt:end),Dpoh(1:dt:end,:),posc,alp,'-',posce);           
            plot(th(1:dt:end),Dprmh(1:dt:end),'--','LineWidth',1,'Color',[pric 1]); % 1
            plot(th(1:dt:end),Dpomh(1:dt:end),'LineWidth',2,'Color',[posc 3*alp]);
        end
        for wyr=1:Ny
            tmp(wyr).FaceAlpha=0.05;
        end
    end

    

    pt(1)=plot(nan,nan,'--','LineWidth',1,'Color',[pric 1]); % 1
    pt(2)=plot(nan,nan,'LineWidth',2,'Color',[posc 3*alp]); % 1


    
    

    

    
    tov=datevec(to); towy=tov(:,1); tom=tov(:,2);
    towy(tom>=9)=towy(tom>=9)+1;
    these=zeros(numel(to),1,'logical');
    
    if any(strcmp(expc,'F'))
        lowys=r.p.wys;
        notthese=lowys>=2000&lowys<=2023;
        lowys=lowys(~notthese);
    else
        lowys=2000:2024;
    end
    for k=1:numel(to)
        if any(lowys==towy(k))
            these(k)=1;
        end
    end

    

    tc=to(~these);
    Dc=Do(~these);
    Do=Do(these);
    to=to(these);

    if expt<3
        falp=1;
        msz=20;%
        lw=0.2;
        ec=[0 0 0];
    else
        falp=1;
        msz=25;
        ec=[0 0 0];
        lw=0.5;
    end

    if expt<3
        plot(t(1:dt:end),Dprm(1:dt:end),'--','LineWidth',1,'Color',[pric 1]); % 0.75 1
        plot(t(1:dt:end),Dpom(1:dt:end),'LineWidth',2,'Color',[posc 3*alp]);
    end
    
    pt(4)=scatter(to,Do,msz,'s','MarkerEdgeColor',ec,...
    'MarkerFaceColor',[0.78 0.41 0.08],'MarkerFaceAlpha',falp,...
    'LineWidth',lw);
    
  
    pt(3)=scatter(tc,Dc,msz,'o','MarkerEdgeColor',ec,...
        'MarkerFaceColor',obsc,'MarkerFaceAlpha',falp,...
        'LineWidth',lw);


   
    if any(strcmp(expc,'F'))
        tvc=datenum('01-Sep-2023');
        plot([tvc tvc],yls,'Color',[0.5 0.5 0.5 0.4],'LineWidth',0.5,'LineStyle','-'); hold on;
    end
    ty=datevec(t); ty=ty(:,1);
    Nticks=(max(ty)-1)-min(ty)+1;
    xticks=zeros(Nticks,1);
    for tks=1:Nticks
        xticks(tks)=datenum(sprintf('01-Sep-%d',min(ty)+tks-1));
    end
    set(gca,'XTick',xticks);

    % Labels are in terms of water year, so WY2000 is the tick at 01.09.1999
    xtls=compose('%d',((min(ty)+1):max(ty))');

    if dolong
        do5=1;
    else
        do5=0;
    end
    if do5
        % Only every 5 years
        xty=datevec(xticks); xty=xty(:,1)+1; % Labels are one year ahead (WY)
        these=~mod(xty,5);
        xty(~these)=missing;
        xtls=compose('%d',xty);
        xtls=strrep(xtls,'NaN','');
    end

    xlim([min(t) max(t)]);
    if ~dolong
        xlim([datenum('01-Sep-1989') datenum('01-Sep-2009')]);
    end
    datetick('x','yyyy','keepticks','keeplimits'); %mm.yy
    axis xy;
    set(gca,'XTick',xticks);
    if j==3
        set(gca,'XTickLabel',xtls);
    else
        set(gca,'XTickLabel',[]);
    end
    xtickangle(0)
    fs=18;
    box on; 
    
    set(groot, 'defaultAxesTickLabelInterpreter','tex');
    set(gca,'TickLabelInterpreter','tex','FontName',fname);
    set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
    tis=sprintf('(%s) %s %s',alphabet(j),fsites{s},methods{j});
    title(tis,'Interpreter','tex','FontSize',1.*fs);

    if j==1
        if expt==3
            if any(strcmp(site,'AB'))
                locis='SouthWest';
            elseif any(strcmp(site,'SG'))
                locis='SouthWest';
            elseif any(strcmp(site,'SB'))
                locis='NorthEast';
            elseif any(strcmp(site,'CD'))
                locis='South';
            else
                locis='NorthEast';
            end
        else
            locis='NorthEast';
        end
        leg=legend(pt(:),{'Prior','Posterior','Calibration','Validation'},...
            'Interpreter','tex','FontSize',0.8*fs,...
            'Location',locis,'Orientation','Horizontal');
    end
    set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.001, 0.0025]);
   
end


linkaxes(findall(gcf, 'type', 'axes'), 'x');
ax=gca;
ax.Toolbar.Visible = 'off';
set(gcf, 'Renderer', 'Painters');

if expt<3
        if ~plotF
            ylabel(tld,'Snow water equivalent [mm w.e.]','Interpreter','tex','FontSize',1.2*fs);
        else
            ylabel(tld,'Fractional snow-covered area [-]','Interpreter','tex','FontSize',1.2*fs)
        end
    else
        ylabel(tld,'Annual cumulative glacier-wide mass balance [mm w.e.]','Interpreter','tex','FontSize',1.2*fs);
    end

pause(2);
if ~plotF
    toprint=sprintf('results/Long%d_Exp_%s_D_%s.pdf',dolong,expc,sites{s});
else
    toprint=sprintf('results/Long%d_Exp_%s_F_%s.pdf',dolong,expc,sites{s});
end

exportgraphics(tld, toprint, 'ContentType', 'vector');

end



    end
end