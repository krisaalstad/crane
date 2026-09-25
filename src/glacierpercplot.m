function flag=glacierpercplot(t,p,obs,Dprm,Dprs,Dpom,Dpos,...
    Cprm,Cprs,Cpom,Cpos,glacier)


% Color settings
pric=[0.8 0 0];
posc=[0 0 0.8];
ec='none';
vc=[0 0 0];
sc=[0.6 0.6 0];
close all;


ty=datevec(t); ty=ty(:,1);
Nticks=(max(ty)-1)-min(ty)+1;
xticks=zeros(Nticks,1);
for j=1:Nticks
    xticks(j)=datenum(sprintf('01-Sep-%d',min(ty)+j-1));
end
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

if p.visD
    fis=figure('units','inch','position',[0,0,18,6]);
else
    fis=figure('units','inch','position',[0,0,18,6]);
end
op=0.05;
nsig=3;

dom=1;
if dom==1
    cf=1e3;
end

Dpri=[Dprm-nsig.*Dprs Dprm Dprm+nsig.*Dprs]./cf;
Dpst=[Dpom-nsig.*Dpos Dpom Dpom+nsig.*Dpos]./cf;


if p.visD

subplot(2,1,1);
pp(1)=percplot(t,Dpri,pric,0.4,'-',ec); hold on;
pp(2)=percplot(t,Dpst,posc,0.4,'-',ec);
pp(1).FaceAlpha=0.1; pp(2).FaceAlpha=0.2;


pt(1)=plot(t,Dprm./cf,'LineWidth',1,'Color',[pric 0.5]);
pt(2)=plot(t,Dpom./cf,'LineWidth',1,'Color',[posc 0.5]);

to=obs.t;
Do=obs.D;

scatter(to,Do./cf,'MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[0.8 0.8 0]);
toy=datevec(to); toy=toy(:,1);
if p.lowy==1
    tov=datevec(to); towy=tov(:,1); tom=tov(:,2);
    towy(tom>=9)=towy(tom>=9)+1;
    these=zeros(numel(to),1,'logical');
    for j=1:numel(to)
        if any(p.lowys==towy(j))
            these(j)=1;
        end
    end
    scatter(to(these),Do(these)./cf,'s','MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[0.8 0.8 0.8]);
end
fs=24;
tis=sprintf('%s',glacier);
title(tis,'Interpreter','Latex','FontSize',fs);
set(gca,'XTick',xticks);

nn=min([min(Do) min(Dpom)]);
xx=max([max(Do) max(Dpom)]);

ylim([1.5.*nn./cf 1.5*xx./cf])
xlim([min(t) max(t)]);
set(gca,'XTick',xticks);
datetick('x','yyyy','keepticks','keeplimits');
set(gca,'XTickLabel',xtls);
xtickangle(0)
axis xy;
box on; 
set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.0015, 0.0025]);
set(groot, 'defaultAxesTickLabelInterpreter','LaTex');
set(gca,'TickLabelInterpreter','LaTex');
set(gca,'GridLineStyle',':','FontSize',0.5.*fs)
ylabel('Annual SMB [m w.e.]','Interpreter','Latex','FontSize',fs);
title(tis,'Interpreter','Latex','FontSize',fs);

end

if p.visD
    subplot(2,1,2);
end
    
Cobsa=[0; cumsum(obs.D(2:2:end),1)];
Cobsw=Cobsa(1:(end-1))+obs.D(1:2:end);
Cobsa=Cobsa(2:end); % Remove initial zero
Cobs=[Cobsw Cobsa]';
Cobs=Cobs(:);

yearinds=(1:(numel(Cobs)/2))';
sigo=zeros(numel(yearinds),2);
sigo(:,2)=sqrt(yearinds.*(2*p.sigD).^2); 
sigo(:,1)=sqrt((yearinds-1).*(2*p.sigD).^2+...
    p.sigD.^2);
% Assume error of 2*p.sigD in AMB and p.sigD in WMB
sigo=sigo'; sigo=sigo(:);

Cop=[Cobs-nsig.*sigo Cobs Cobs+nsig.*sigo]./cf;
Cop1=[Cobs-sigo Cobs Cobs+sigo]./cf;


Cpri=[Cprm-nsig.*Cprs Cprm Cprm+nsig.*Cprs]./cf;
Cpri1=[Cprm-Cprs Cprm Cprm+Cprs]./cf;
Cpst=[Cpom-nsig.*Cpos Cpom Cpom+nsig.*Cpos]./cf;
Cpst1=[Cpom-Cpos Cpom Cpom+Cpos]./cf;

percplot(t,Cpri,pric,0.1,'-',ec); hold on;
pp(1)=percplot(t,Cpri1,pric,0.1,'-',ec); hold on;
percplot(t,Cpst,posc,0.1,'-',ec);
pp(2)=percplot(t,Cpst1,posc,0.1,'-',ec);

percplot(obs.t,Cop,[0.2 0.2 0],0.1,'-',ec);
percplot(obs.t,Cop1,[0.2 0.2 0],0.1,'-',ec);

pt(1)=plot(t,Cprm./cf,'LineWidth',1,'Color',[pric 0.5]);
pt(2)=plot(t,Cpom./cf,'LineWidth',1,'Color',[posc 0.5]);


plot(obs.t,Cobs./cf,'Color',[0 0 0]);
if p.lowy==1
    to=obs.t;
    tov=datevec(to); towy=tov(:,1); tom=tov(:,2);
    towy(tom>=9)=towy(tom>=9)+1;
    these=zeros(numel(to),1,'logical');
    for j=1:numel(to)
        if any(p.lowys==towy(j))
            these(j)=1;
        end
    end
    pt(4)=scatter(to(these),Cobs(these)./cf,'s','MarkerEdgeColor',[0 0 0],...
    'MarkerFaceColor',[0.8 0.8 0.8]);
    pt(3)=scatter(obs.t(~these),Cobs(~these)./cf,'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',[0.8 0.8 0]);
else

    plot(obs.t,Cobs./cf,'Color',[0 0 0]);
    scatter(obs.t,Cobs./cf,'MarkerEdgeColor',[0 0 0],...
        'MarkerFaceColor',[0.8 0.8 0]);
end


fs=24;
ylabel('Cumulative SMB [m w.e.]','Interpreter','Latex','FontSize',fs);
set(gca,'XTick',xticks);
xlim([min(t) max(t)]);
datetick('x','yyyy','keepticks','keeplimits'); %mm.yy
set(gca,'XTickLabel',xtls);
xtickangle(0)
axis xy;
box on; 
set(gca,'TickDir','out','LineWidth',1.5,'TickLength',[0.0015, 0.0025]);
set(groot, 'defaultAxesTickLabelInterpreter','LaTex');
set(gca,'TickLabelInterpreter','LaTex');
set(gca,'GridLineStyle',':','FontSize',0.5.*fs);
ylabel('Cumulative MB [m w.e.]','Interpreter','Latex','FontSize',fs);

alphabet = 'abcdefghijklmonpqrstuvwxyz';
if any(strcmp(p.case,'CP'))
    leg=legend(pt(:),{'Prior','Posterior','Calibration','Validation'},...
        'Interpreter','Latex','FontSize',fs,...
        'Location','SouthWest');
end
if any(strcmp(glacier,'AB'))
    lstart=1;
    yls=[-50 5];
elseif any(strcmp(glacier,'SG'))
    lstart=4;
    yls=[-50 5];
elseif any(strcmp(glacier,'SB'))
    lstart=7;
    yls=[-50 5];
elseif any(strcmp(glacier,'CD'))
    lstart=10;
    yls=[-45 10];
end
letter=alphabet(lstart);
if any(strcmp(p.case,'NP'))
    letter=alphabet(lstart+1);
elseif any(strcmp(p.case,'PP'))
    letter=alphabet(lstart+2);
elseif any(strcmp(p.case,'MAP'))
    letter=alphabet(lstart+3);
end

if ~p.visD
    tis=sprintf('%s %s',glacier,p.case);
    title(tis,'Interpreter','Latex','FontSize',fs);
end


nn=min([min(Cobs) min(Cpom)]);
xx=max([max(Cobs) max(Cpom)]);

ylim(yls)

ax=gca;
ax.Toolbar.Visible = 'off';


pause(2);
toprint=sprintf('results/Glacier_%s_LO%d_%s.jpg',glacier,...
    p.lowy,p.case);
exportgraphics(gcf,toprint,'Resolution',300);

flag=1;

end