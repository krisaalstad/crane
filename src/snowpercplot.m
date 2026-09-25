function flag=snowpercplot(t,p,obs,Dprm,Dprs,Dpom,Dpos,...
    Fprm,Fprs,Fpom,Fpos,snowsite)


% Color settings
pric=[0.8 0 0];
posc=[0 0 0.8];
ec='none';
mc=[0.5 0.5 0.5];%'none';
vc=[0 0 0];
sc=[0.6 0.6 0];
close all;


if p.assimF
    fis=figure('units','inch','position',[0,0,18,12]);
else
    fis=figure('units','inch','position',[0,0,18,6]);
end


op=0.05;
nsig=1;

dom=1;
if dom==1
    cf=1e3;
end

Dpri=[Dprm-nsig.*Dprs Dprm Dprm+nsig.*Dprs]./cf;
Dpst=[Dpom-nsig.*Dpos Dpom Dpom+nsig.*Dpos]./cf;
Dpri=max(Dpri,0);
Dpst=max(Dpst,0);

if ~p.vispri
    Dpri=zeros(size(Dpri));
    Dprm=zeros(size(Dprm));
end

if ~isempty(Fpom)
    subplot(2,1,1);
end
pp(1)=percplot(t,Dpri,pric,0.1,'-',ec); hold on;
pp(2)=percplot(t,Dpst,posc,0.1,'-',ec);
pp(1).FaceAlpha=0.1; pp(2).FaceAlpha=0.1;
to=obs.t;
Do=obs.D;

if ~p.assimD
    Doc=[0.8 0.8 0.8];
else
    Doc=[0.8 0.8 0];
end


if p.assimD&&p.lowy
    tov=datevec(to); towy=tov(:,1); tom=tov(:,2);
    towy(tom>=9)=towy(tom>=9)+1;
    these=zeros(numel(to),1,'logical');
    for j=1:numel(to)
        if any(p.lowys==towy(j))
            these(j)=1;
        end
    end
    tc=to(~these);
    Dc=Do(~these);
    Do=Do(these);
    to=to(these);
    scatter(tc,Dc./cf,20,'o','MarkerEdgeColor',mc,...
    'MarkerFaceColor',[0.8 0.8 0],'MarkerFaceAlpha',0.3);
end

% Dummy plots:
pt(1)=scatter([],[],20,'s','MarkerEdgeColor',mc,...
    'MarkerFaceColor',[0.8 0.8 0.8],'MarkerFaceAlpha',0.3);
pt(3)=scatter([],[],20,'o','MarkerEdgeColor',mc,...
    'MarkerFaceColor',[0.8 0.8 0],'MarkerFaceAlpha',0.3);


pt(4)=scatter(to,Do./cf,20,'s','MarkerEdgeColor',mc,...
    'MarkerFaceColor',[0.8 0.8 0.8],'MarkerFaceAlpha',0.3);

pt(1)=plot(t,Dprm./cf,'LineWidth',1,'Color',[pric 0.6]);
pt(2)=plot(t,Dpom./cf,'LineWidth',1,'Color',[posc 0.6]);

fs=24;

ty=datevec(t); ty=ty(:,1);
Nticks=(max(ty)-1)-min(ty)+1;
xticks=zeros(Nticks,1);
for j=1:Nticks
    xticks(j)=datenum(sprintf('01-Sep-%d',min(ty)+j-1));
end
set(gca,'XTick',xticks);
% Labels are in terms of water year, so WY2000 is the tick at 01.09.1999
xtls=compose('%d',((min(ty)+1):max(ty))');

do5=1;
if do5
    % Only every 5 years
    xty=datevec(xticks); xty=xty(:,1)+1; % Labels are one year ahead (WY)
    these=~mod(xty,5);
    xty(~these)=missing;
    xtls=compose('%d',xty);
    xtls=strrep(xtls,'NaN','');
    
end



xlim([min(t) max(t)]);
datetick('x','yyyy','keepticks','keeplimits'); %mm.yy
axis xy;
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xtls);
xtickangle(0)
box on; 
set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.0015, 0.0025]);
set(groot, 'defaultAxesTickLabelInterpreter','LaTex');
set(gca,'TickLabelInterpreter','LaTex');
set(gca,'GridLineStyle',':','FontSize',0.5.*fs)



alphabet = 'abcdefghijklmonpqrstuvwxyz';
if any(strcmp(snowsite,'HS'))
    lstart=1;
    if any(strcmp(p.case,'CP'))
        leg=legend(pt(:),{'Prior','Posterior','Assimilation','Validation'},...
            'Interpreter','Latex','FontSize',fs,...
            'Location','NorthEast','Orientation','Horizontal');
    end
    yls=[0 0.5];
elseif any(strcmp(snowsite,'AK'))
    lstart=4;
    if any(strcmp(p.case,'CP'))
        leg=legend(pt(:),{'Prior','Posterior','Assimilation','Validation'},...
            'Interpreter','Latex','FontSize',fs,...
            'Location','NorthEast','Orientation','Horizontal');
    end
    yls=[0 0.5];
elseif any(strcmp(snowsite,'FF'))
    lstart=7;
    if any(strcmp(p.case,'CP'))
        leg=legend(pt(:),{'Prior','Posterior','Assimilation','Validation'},...
            'Interpreter','Latex','FontSize',fs,...
            'Location','NorthEast','Orientation','Horizontal');
    end
    yls=[0 0.7];
elseif any(strcmp(snowsite,'WF'))
    lstart=10;
    if any(strcmp(p.case,'CP'))
        leg=legend(pt(:),{'Prior','Posterior','Assimilation','Validation'},...
            'Interpreter','Latex','FontSize',fs,...
            'Location','NorthEast','Orientation','Horizontal');
    end
    yls=[0 1.4];
end
letter=alphabet(lstart);
if any(strcmp(p.case,'NP'))
    letter=alphabet(lstart+1);
elseif any(strcmp(p.case,'PP'))
    letter=alphabet(lstart+2);
elseif any(strcmp(p.case,'MAP'))
    letter=alphabet(lstart+3);
end
ylim(yls)
tis=sprintf('%s %s',snowsite,p.case);
title(tis,'Interpreter','Latex','FontSize',fs);
ylabel('SWE [m w.e.]','Interpreter','Latex','FontSize',fs);

nsig=1;
Fpri=[Fprm-nsig.*Fprs Fprm Fprm+nsig.*Fprs];

Fpri=max(Fpri,0);
Fpri=min(Fpri,1);

Fpst=[Fpom-nsig.*Fpos Fpom Fpom+nsig.*Fpos];
Fpst=max(Fpst,0);
Fpst=min(Fpst,1);

if ~isempty(Fpom)

tos=obs.ts;
Fos=obs.Fs;
    if ~p.vispri
        Fpri=zeros(size(Fpri));
        Fprm=zeros(size(Fprm));
    end

    subplot(2,1,2);
pp(1)=percplot(t,Fpri,pric,0.1,'-',ec); hold on;
pp(2)=percplot(t,Fpst,posc,0.1,'-',ec);
pp(1).FaceAlpha=0.1; pp(2).FaceAlpha=0.2;

pt(1)=plot(t,Fprm,'LineWidth',1,'Color',[pric 0.5]);
pt(2)=plot(t,Fpom,'LineWidth',1,'Color',[posc 0.5]);


scatter(tos,Fos,20,'MarkerEdgeColor',mc,...
    'MarkerFaceColor',[0.8 0.8 0]);%,'MarkerFaceAlpha',0.2);
toy=datevec(tos); toy=toy(:,1);
if p.lowy==1
    tov=datevec(tos); towy=tov(:,1); tom=tov(:,2);
    towy(tom>=9)=towy(tom>=9)+1;
    these=zeros(numel(tos),1,'logical');
    for j=1:numel(tos)
        if any(p.lowys==towy(j))
            these(j)=1;
        end
    end
    scatter(tos(these),Fos(these),'MarkerEdgeColor',mc,...
    'MarkerFaceColor',[0.8 0.8 0.8]);
end

fs=24;
ylabel('fSCA [-]','Interpreter','Latex','FontSize',fs);

Nticks=(2021-1999+1); % Location of ticks
xticks=zeros(Nticks,1);
for j=1:Nticks
    xticks(j)=datenum(sprintf('01-Sep-%d',1998+j));
end
set(gca,'XTick',xticks);
% Labels are in terms of water year, so WY2000 is the tick at 01.09.1999
xtls=compose('%d',(2000:2022)');
tses=datenum('01-Sep-1999'); % Start of MODIS era
tsee=datenum('01-Sep-2022'); 



xlim([tses tsee]);
datetick('x','yyyy','keepticks','keeplimits'); 
axis xy;
set(gca,'XTickLabel',xtls);
box on; 
set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.0015, 0.0025]);
set(groot, 'defaultAxesTickLabelInterpreter','LaTex');
set(gca,'TickLabelInterpreter','LaTex');
set(gca,'GridLineStyle',':','FontSize',0.5.*fs);

ylabel('fSCA [-]','Interpreter','Latex','FontSize',fs);
ylim([0 1])

end

ax=gca;
ax.Toolbar.Visible = 'off';


pause(2);
toprint=sprintf('results/Snow_%s_D%d_F%d_LO%d_%s.jpg',snowsite,...
    p.assimD,p.assimF,p.lowy,p.case);
exportgraphics(gcf,toprint,'Resolution',300);


flag=1;

end
