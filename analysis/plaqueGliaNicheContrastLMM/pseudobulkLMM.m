%% Check replication with subject
clear all;
close all;

addpath(genpath('gseaMATLAB/'));

%% Parameters
stratA = 'low_RGN_pos'; % see header of pseudobulk.tsv
stratB = 'low_RGN_neg';
layerSel = 'GM';

%% Load pseudobulk expression
fid = fopen('pseudobulk.tsv');
header = regexp(fgetl(fid),'\t','split');
temp = textscan(fid,['%s',repmat('%f',1,length(header)-1)],'Delimiter','\t');
fclose(fid);
subj = header(2:end)';
gene = temp{1};
nGene = length(gene);
feat = cell2mat(temp(2:end))';
        
%% Split by condition
indA = contains(subj,stratA) & contains(subj,layerSel);
subjA = subj(indA);
featA = feat(indA,:);
indB = contains(subj,stratB) & contains(subj,layerSel);
subjB = subj(indB);
featB = feat(indB,:);

%% Concatenate A and B
subj = [subjA;subjB];
label = [ones(length(subjA),1);zeros(length(subjB),1)];
feat = [featA;featB];
subj = regexprep(subj,'_(\w+)',''); % removes annotations 
nSect = length(subj);
        
%% Load disease diagnosis
fid = fopen('diagnosis.csv');
header = regexp(fgetl(fid),',','split');
temp = textscan(fid,'%s%f%f','Delimiter',',');
fclose(fid);
[~,ind] = ismember(subj,temp{1});
ad = temp{3}(ind);
        
%% Confounds
fid = fopen('confounds.csv');
header = regexp(fgetl(fid),',','split');
temp = textscan(fid,'%s%f%f%f','Delimiter',',');
fclose(fid);
[~,ind] = ismember(subj,temp{1});
rin = temp{2}(ind);
libBatch = temp{3}(ind);
age = temp{4}(ind);
        
%% Convert library Batch to binary
[c,~,indC] = unique(libBatch);
nC = length(c)-1;
lb = zeros(nSect,nC);
for l = 1:nC
    lb(indC==l,l) = 1;
end
indStd = std(lb)>0;
lb = lb(:,indStd);
        
%% Statistical test
subj = categorical(subj);
beta = nan(nGene,1);
se = nan(nGene,1);
t = nan(nGene,1);
p = nan(nGene,1);
parfor k = 1:nGene
    if mean(ad)<1
        X = [ones(nSect,1),age,rin,lb,ad,label];
        if sum(~isnan(feat(:,k)))>2*size(X,2) && rank(X)==size(X,2)
            lme = fitlmematrix(X,feat(:,k),{ones(nSect,1)},{subj});
            [~,~,stats] = fixedEffects(lme);
            beta(k) = stats.Estimate(end);
            se(k) = stats.SE(end);
            t(k) = stats.tStat(end);
            p(k) = stats.pValue(end);
        end
    else % Some conditions might have only AD subjects
        X = [ones(nSect,1),age,rin,lb,label];
        if sum(~isnan(feat(:,k)))>2*size(X,2) && rank(X)==size(X,2)
            lme = fitlmematrix(X,feat(:,k),{ones(nSect,1)},{subj});
            [~,~,stats] = fixedEffects(lme);
            beta(k) = stats.Estimate(end);
            se(k) = stats.SE(end);
            t(k) = stats.tStat(end);
            p(k) = stats.pValue(end);
        end
    end
end

%% GO geneset
[gsGO,pGO,nesGO] = gsea(t,gene,'gseaMATLAB/c5.all.v2022.1.Hs.symbols.gmt');
[~,indS] = sort(pGO);
table(gsGO(indS(1:10)),pGO(indS(1:10))) % Exemplar GO genesets

%% Active glia geneset
[gs,p,nes] = gsea(t,gene,'gseaMATLAB/activeGlia.gmt');
table(gs,p)

%% Sun et al. MG geneset
[gs,p,nes] = gsea(t,gene,'gseaMATLAB/mgSun2023.gmt');
table(gs,p)

%% RUSH Abeta response
[gs,p,nes] = gsea(t,gene,'gseaMATLAB/rushSCabResponse.gmt');
table(gs,p)

%% RUSH MG Abeta response
[gs,p,nes] = gsea(t,gene,'gseaMATLAB/rushSCmgAbResponse.gmt');
table(gs,p)


