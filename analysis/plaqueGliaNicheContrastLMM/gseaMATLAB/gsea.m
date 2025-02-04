%% GSEA
% Input:    t = t-values of genes
%           gene = list of genes
%           gsPath = geneset path with geneset filename
%           std = 0 or 1
% Output:   geneset = tested genesets
%           p = p-value of geneset enrichment
%           nes = normalized enrichment score
function [geneset,p,nes] = gsea(t,gene,gsPath,std)
fnPath = 'gseaMATLAB/';
if nargin < 4
    std = 1;
end

%% Split geneset path
[gsPath,gsFile,gsExt] = fileparts(gsPath);

%% Generate ranked gene scores
[~,indS] = sort(t,1,'descend');
pS = t(indS);
geneS = gene(indS);
T = table(geneS,pS);
rn = randi(1e6,1);
tfolder = ['temp',int2str(rn),'/'];
while exist(tfolder,'dir')
    rn = randi(1e6,1);
    tfolder = ['temp',int2str(rn),'/'];
end
system(['mkdir ',tfolder]);
writetable(T,[tfolder,'rankedPval.txt'],'Delimiter','\t'); 
system(['mv ',tfolder,'rankedPval.txt ',tfolder,'rankedPval.rnk']);
system(['Rscript ',fnPath,'gsea.R ',gsPath,'/ ',gsFile,gsExt,' ',tfolder]);
pause(5);
system(['rm ',tfolder,'rankedPval.rnk']);

%% GSEA
if std==1
    fid = fopen([tfolder,gsFile,'.gsea.std.txt']);
    header = regexp(fgetl(fid),'\t','split');
    temp = textscan(fid,'%s%f%f%f%f%f%f%s','Delimiter','\t');
    fclose(fid);
    geneset = temp{1};
    p = temp{2};
    nes = temp{6};
else
    fid = fopen([tfolder,gsFile,'.gsea.pos.txt']);
    header = regexp(fgetl(fid),'\t','split');
    temp = textscan(fid,'%s%f%f%f%f%f%f%s','Delimiter','\t');
    fclose(fid);
    geneset = temp{1};
    pPos = temp{2};
    nesPos = temp{6};
    fid = fopen([tfolder,gsFile,'.gsea.neg.txt']);
    header = regexp(fgetl(fid),'\t','split');
    temp = textscan(fid,'%s%f%f%f%f%f%f%s','Delimiter','\t');
    fclose(fid);
    geneset = temp{1};
    pNeg = temp{2};
    nesNeg = temp{6};
    nSet = length(geneset);
    p = nan(nSet,1);
    nes = nan(nSet,1);
    for i = 1:nSet
        [p(i),indMin] = min([pPos(i),pNeg(i)]);
        if indMin==1
            nes(i) = nesPos(i);
        else
            nes(i) = nesNeg(i);
        end
    end
    p = min(2*p,1); % 2 tail test
    nes = nes/1.5; % Account for 2 tails, 1.5 is empirical
end
% system(['rm ',gsFile,'.gsea.*.txt']);
system(['rm -r ',tfolder]);
