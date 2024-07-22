initCobraToolbox
load('recon1.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Not lethal in section 5 and lethal in section 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do single gene deletion for all on the model that gene is not lethal -
%cholangiocytes
df = readcell('PM_genes_indicator_for_singleGeneDeletion.csv');
dfs = string(df); %transform data type to string for easier indexing


dfunique = unique(dfs, 'stable'); %get unique genes

%WT growth rate for section 5
load('mouse_intestine_relative_section5_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat');
optimizeCbModel(PM);
grRateWT = ans.f;

outchol = cell(length(dfunique),1);
for x = 1:length(dfunique)

[temp,lethal_hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(PM,PM.genes(dfs == dfunique(x)));
optimizeCbModel(temp); 
grRateKO = ans.f;
grRatio = grRateKO/grRateWT;

outchol{x,:} = grRatio;
end

df = readtable('int_OptimizeCbModel_diff_lethal_in_relative_section1.csv');
df1 = string(df.NCBI); %transform data type to string for easier indexing
df_hgnc = df.hgnc;

for x = 1:length(df1)
%%gene knockout lethal cell type
load('mouse_intestine_relative_section1_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat');

[temp,lethal_hasEfffect,hepconstrRxnNames,deletedGenes] = deleteModelGenes(PM,PM.genes(dfs == df1(x)));
Hep_lethal_biomass = [biomassPrecursorCheck(temp)];

%%gene knockout not lethal cell type
%load metabolic model 
load('mouse_intestine_relative_section5_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat');
%load list of genes that are not lethal here but lethal in other cell type

%delete gene of interest
[new,hasEffect,cholconstrRxnNames,deletedGenes] = deleteModelGenes(PM,PM.genes(dfs == df1(x)));
chol_biomass_in_hep_lethal = [biomassPrecursorCheck(new)];

%remove related reactions in model %use PM
temp=removeRxns(PM,cholconstrRxnNames);


%get growth rate ratio of secondary knock out

    outchol2 = cell(length(dfunique),1);
    for y = 1:length(dfunique)
    %remove related reactions in model %use PM
    temp=removeRxns(PM,cholconstrRxnNames);
    [temp,lethal_hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(temp,temp.genes(dfs == dfunique(y)));%try singleRxnDeletion
    optimizeCbModel(temp); 
    grRateKO2 = ans.f;
    grRatio = grRateKO2/grRateWT;

    outchol2{y,:} = grRatio;
    end
     
%alternative2be=dfunique(string(outchol)>'0.8' & string(outchol2)<'0.2');
%let's do this in R
           
out = [cellstr(dfunique) outchol outchol2];%for lethal genes in secondary KO, find biomass?
           
%alt_res = [hasEffect;alternative2be];

str = '(lethal_in_sec1)sec1_lethal_constRxnNames_for_sec1KO.txt';
writecell(hepconstrRxnNames, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec1)sec1_lethal_constRxnNames_for_sec5.txt';
writecell(cholconstrRxnNames, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec1)sec1_lethal_has_effect_for_sec1KO.txt';
writematrix(lethal_hasEfffect, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec1)sec1_lethal_biomass_for_sec1.txt';
writecell(Hep_lethal_biomass, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec1)sec1_lethal_biomass_for_sec5.txt';
writecell(chol_biomass_in_hep_lethal, insertAfter(str,'for_',df_hgnc(x)))

%str = '(lethal_in_hep)alternative_genes_for_chol.txt';
%writematrix(alt_res, insertAfter(str,'for_',df_hgnc(x)));
str = '(lethal_in_sec1)grRatio_KO_for_sec5';
writecell(out, insertAfter(str,'grRatio_',df_hgnc(x)));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Not lethal in section 1 and lethal in section 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%do single gene deletion for all on the model that gene is not lethal -
%hep 1
df = readcell('PM_genes_indicator_for_singleGeneDeletion.csv');
dfs = string(df); %transform data type to string for easier indexing


dfunique = unique(dfs, 'stable'); %get unique genes

%WT growth rate for section 1
load('mouse_intestine_relative_section1_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat');
optimizeCbModel(PM);
grRateWT = ans.f;

outhep = cell(length(dfunique),1);
for x = 1:length(dfunique)

[temp,lethal_hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(PM,PM.genes(dfs == dfunique(x)));
optimizeCbModel(temp); 
grRateKO = ans.f;
grRatio = grRateKO/grRateWT;

outhep{x,:} = grRatio;
end

df = readtable('int_OptimizeCbModel_diff_lethal_in_relative_section5.csv');
df1 = string(df.NCBI); %transform data type to string for easier indexing
df_hgnc = df.hgnc;

for x = 1:length(df1)
%%gene knockout lethal cell type
load('mouse_intestine_relative_section5_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat');

[temp,lethal_hasEfffect,cholconstrRxnNames,deletedGenes] = deleteModelGenes(PM,PM.genes(dfs == df1(x)));
Chol_lethal_biomass = [biomassPrecursorCheck(temp)];

%%gene knockout not lethal cell type
%load metabolic model 
load('mouse_intestine_relative_section1_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat');
%load list of genes that are not lethal here but lethal in other cell type

%delete gene of interest
[new,hasEffect,hepconstrRxnNames,deletedGenes] = deleteModelGenes(PM,PM.genes(dfs == df1(x)));
hep_biomass_in_chol_lethal = [biomassPrecursorCheck(new)];

%remove related reactions in model %use PM
temp=removeRxns(PM,hepconstrRxnNames);


%get growth rate ratio of secondary knock out

    outhep2 = cell(length(dfunique),1);
    for y = 1:length(dfunique)
    
    %remove related reactions in model %use PM
    temp=removeRxns(PM,hepconstrRxnNames);
    [temp,lethal_hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(temp,temp.genes(dfs == dfunique(y)));
    optimizeCbModel(temp); 
    grRateKO2 = ans.f;
    grRatio = grRateKO2/grRateWT;

    outhep2{y,:} = grRatio;
    end
     
%alternative2be=dfunique(string(outhep)>'0.8' & string(outhep2)<'0.2');

           
out = [cellstr(dfunique) outhep outhep2];
           
%alt_res = [hasEffect;alternative2be];

str = '(lethal_in_sec5)sec5_lethal_constRxnNames_for_sec5KO.txt';
writecell(cholconstrRxnNames, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec5)sec5_lethal_constRxnNames_for_sec1.txt';
writecell(hepconstrRxnNames, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec5)sec5_lethal_has_effect_for_sec5KO.txt';
writematrix(lethal_hasEfffect, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec5)sec5_lethal_biomass_for_sec5.txt';
writecell(Chol_lethal_biomass, insertAfter(str,'for_',df_hgnc(x)))
str = '(lethal_in_sec5)sec5_lethal_biomass_for_sec1.txt';
writecell(hep_biomass_in_chol_lethal, insertAfter(str,'for_',df_hgnc(x)))

%str = '(lethal_in_chol)alternative_genes_for_hep.txt';
%writematrix(alt_res, insertAfter(str,'for_',df_hgnc(x)));
str = '(lethal_in_sec5)grRatio_KO_for_sec1';
writecell(out, insertAfter(str,'grRatio_',df_hgnc(x)));
end