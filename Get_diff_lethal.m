%read cell without gene variant information
df = readcell('PM_genes_indicator_for_singleGeneDeletion.csv');
dfs = string(df); %transform data type to string for easier indexing

%WT growth rate for section1
load('mouse_intestine_relative_section1_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat')
optimizeCbModel(PM);
grRateWT = ans.f;

dfunique = unique(dfs, 'stable'); %get unique genes

outhep = cell(length(dfunique),1);
for x = 1:length(dfunique)

[temp,lethal_hasEffect,constrRxnNames,deletedGenes] = deleteModelGenes(PM,PM.genes(dfs == dfunique(x)));
optimizeCbModel(temp); 
grRateKO = ans.f;
grRatio = grRateKO/grRateWT;

outhep{x,:} = grRatio;
end

%WT growth rate for cholangiocytes
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

hepchol = [cellstr(dfunique) outhep outchol];

str = 'Int_deleteModelGenes_optimizeCbModel_for_relative_section1and5.txt';
writecell(hepchol, str)

%data is read in R for further processing to get differentially lethal
%genes for each cell type