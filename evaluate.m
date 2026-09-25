clearvars;
close all;
tardir='results/';

% Check input data exists, otherwise fetch it from zenodo
if ~any(exist(tardir,'dir'))
    get_results;
end

set(0, 'DefaultFigureRenderer', 'painters');


ssitestrs={'Hornsund';'Abisko';'Filefjell';'Weissfluhjoch'};
gsitestrs={'A. Brøggerbreen';'Storglaciären';'Storbrean';'Claridenfirn'};
methods={'Prior'; 'CP'; 'NP'; 'PP'; 'MAP'; 'PMCMC'};
ssites={'HS'; 'AK'; 'FF'; 'WF'};
gsites={'AB'; 'SG'; 'SB'; 'CD'};
experiments={'D1_F0_LO1'; 'D0_F1_LO0'; 'LO1'};
compare={'sp','cp','np','hp','pp','map','pmcmc'};
Nexp=3;
Nsites=4;
crpsv=zeros(Nsites*Nexp,numel(methods));
crpsc=zeros(size(crpsv));

%% Loop for computing evaluation metrics
for ex=1:3
    expis=experiments{ex};
    for s=1:Nsites
        fprintf('\n Exp=%d, site=%d \n',ex,s);
        rowis=(ex-1)*Nsites+s;
        if ex<3
            site=ssites{s};
            tarf=sprintf('%s/Snow_%s_%s.mat',...
                tardir,site,expis);
            sitestrs=ssitestrs;
            
        else
            site=gsites{s};
            tarf=sprintf('%s/Glacier_%s_%s.mat',...
                tardir,site,expis);
            sitestrs=gsitestrs;
        end
        load(tarf);



        % Add line to read original obs for SWE experiment (not thinned)
        if ex==1||ex==2
            obsfile=sprintf('input/snow/%s/%s_obs.mat',site,site);
            load(obsfile);
            disp('assigning obs')
            r.obs.t=obs.t; % Overwrite with non thinned obs
            r.obs.D=obs.D;
        end
        

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

        if any(strcmp(expis,'D0_F1_LO0'))
            theseoa=wyo>=2000&wyo<=2023;
        else
            theseoa=wyo<r.p.lowys(1);
        end

        x=r.obs.D;
        for m=1:numel(methods)
            if m==1
                ym=r.pp.Dprim(theseo);  
                ys=r.pp.Dpris(theseo);
            elseif m==2
                ym=r.cp.Dm(theseo); 
                ys=r.cp.Ds(theseo);
            elseif m==3
                ym=r.np.Dm(theseo); 
                ys=r.np.Ds(theseo);
            elseif m==4
                ym=r.pp.Dm(theseo); 
                ys=r.pp.Ds(theseo);
            elseif m==5
                ym=r.map.Dm(theseo); 
                ys=r.map.Ds(theseo);
            elseif m==6
                ym=r.mcmc.Dm(theseo);
                ys=r.mcmc.Ds(theseo);
            end

            expsite(rowis).method(m).crpsc=CRPSg(ym(theseoa),ys(theseoa),x(theseoa));
            expsite(rowis).method(m).crpsv=CRPSg(ym(~theseoa),ys(~theseoa),x(~theseoa));
            crpsc(rowis,m)=mean(CRPSg(ym(theseoa),ys(theseoa),x(theseoa)));
            crpsv(rowis,m)=mean(CRPSg(ym(~theseoa),ys(~theseoa),x(~theseoa)));
           


        end




    end
end

nmt=4; % Number of methods to compare in table (exclude benchmarks)
thesem=1:nmt;
dofi=0;


%% CRPS table
fprintf('\n \n CRPS TABLE \n \n')
for j=1:3
    Ns=4; % Number of sites per experiment

    % Header
    if j==1
        expname=' In situ SWE';
        sites=ssites;
        sitestrs=ssitestrs;
    elseif j==2
        expname=' MODIS FSCA';
        sites=ssites;
        sitestrs=ssitestrs;
    else
        expname=' Mass balance';
        sites=gsites;
        sitestrs=gsitestrs;
    end
    stris=['\\hline' expname];
    for i=1:numel(methods)
        stris=[stris ' & ' methods{i}];
    end
    stris=[stris '\\\\ \\hline \n'];
    fprintf(stris);
    
    FI=zeros(6,2);

    for k=1:numel(sites)
        stris=[sitestrs{k}];
        kis=(j-1)*Ns+k;
        inds=1:nmt;
        indsminc=inds(round(crpsc(kis,thesem))==min(round(crpsc(kis,thesem))));
        indsminv=inds(round(crpsv(kis,thesem))==min(round(crpsv(kis,thesem))));

        

        for i=1:numel(methods)
            if any(indsminv==i)
                textv='\\underline';
            else
                textv='\\mathrm';
            end
            if any(indsminc==i)
                textc='\\mathbf';
            else
                textc='\\mathrm';
            end
            stris=[stris ' & ' sprintf('$%s{%d}(%s{%d})$',...
                textc,round(crpsc(kis,i)),textv,round(crpsv(kis,i)))];
            FI(i,1)=FI(i,1)+crpsc(kis,i)/crpsc(kis,1);
            FI(i,2)=FI(i,2)+crpsv(kis,i)/crpsv(kis,1);
        end
        stris=[stris ' \\\\ \n'];

        strfi=[];
        
        % Printing the FI string
        if k==numel(sites)
            strfi='Mean PI (\\%%)';
            FI=1-FI./numel(sites);
            inds=1:nmt;
            indsmaxc=inds(FI(1:4,1)==max(FI(1:4,1)));
            indsmaxv=inds(FI(1:4,2)==max(FI(1:4,2)));
            for i=1:numel(methods)
                if any(indsmaxv==i)
                    textv='\\underline';
                else
                    textv='\\mathrm';
                end
                if any(indsmaxc==i)
                    textc='\\mathbf';
                else
                    textc='\\mathrm';
                end
                strfi=[strfi ' & ' sprintf('$%s{%d}(%s{%d})$',...
                textc,round(1e2*FI(i,1)),textv,round(1e2*FI(i,2)))];
            end
            strfi=[strfi ' \\\\ \n'];

        end

        
        fprintf(stris);
        if k==numel(sites)
            if j==3
                strfi=[strfi '\\hline'];
            end
            fprintf(strfi);
        end
    end

end