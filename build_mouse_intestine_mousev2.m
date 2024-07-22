%changed line 18 in map_gene_scores_to_rxns, so that reaction with 0
%expression evidence are not penalized - because when using only 1 spot, 0
%expression can also be due to drop out, not true 0
%This is a model for a spot in the PIN region defined by high NPY
%expression.

%current cutoff for "section_1_2_model_mean_new_conn.mat" is 0.2
%fixed a bug, sum([10 NaN]) will be nan, which is then changed to 0. used
%nansum instead. this changed number of core reactions a lot


initCobraToolbox
%load('mouseModel.mat')
load('MouseModel_New.mat')
model = MouseModel_New;
model = update_mouse_model(model);

model=addReaction(model,'biomassRxn','0.60552 ala-L[c] + 0.030122 amp[c] + 0.39141 arg-L[c] + 0.33306 asn-L[c] + 0.26882 asp-L[c] + 30 atp[c] + 0.026219 chsterol[c] + 0.0019665 clpn_hs[c] + 0.030122 cmp[c] + 0.02 crn[c] + 0.16301 cys-L[c] + 0.003769 dag_hs[c] + 0.013565 damp[c] + 0.009074 dcmp[c] + 0.009074 dgmp[c] + 0.013565 dtmp[c] + 0.32832 gln-L[c] + 0.39909 glu-L[c] + 0.51919 gly[c] + 0.34069 glygn2[c] + 0.050294 gmp[c] + 0.01 gthrd[c] + 30 h2o[c] + 0.0074453 hdca[c] + 0.0029781 hdcea[c] + 0.14219 his-L[c] + 0.32929 ile-L[c] + 0.56489 leu-L[c] + 0.0032774 lpchol_hs[c] + 0.57111 lys-L[c] + 0.003769 mag_hs[c] + 0.13315 met-L[c] + 0.003 nad[c] + 0.0003 nadp[c] + 0.0044672 ocdca[c] + 0.025224 ocdcea[c] + 0.030152 pa_hs[c] + 0.0098323 pail_hs[c] + 0.0025 paps[c] + 0.047195 pchol_hs[c] + 0.044245 pe_hs[c] + 0.17578 phe-L[c] + 0.24458 pro-L[c] + 0.015076 ps_hs[c] + 0.0014 ptrc[c] + 4.8648e-05 q10[m] + 0.3973 ser-L[c] + 0.0091768 sphmyln_hs[c] + 0.016 spmd[c] + 0.006 sprm[c] + 7e-05 thbpt[c] + 0.38916 thr-L[c] + 0.041207 trp-L[c] + 0.15451 tyr-L[c] + 0.056957 ump[c] + 0.39711 val-L[c] + 0.022286 xolest_hs[c]  -> 30 adp[c] + 30 h[c] + 30 pi[c]');
model=addReaction(model,'IPDPxc','ipdp[x] 	->	ipdp[c]');
model=addReaction(model,'TYRtm','tyr-L[c] 	->	tyr-L[m]');
model=addReaction(model,'Tyr-ggnt','Tyr-ggn[e] 	->	Tyr-ggn[c]');
model=addReaction(model,'lystopeplys','lys-L[e] 	->	peplys[e]');
model=addReaction(model,'pepslys_sink','pepslys[r] 	->');
ex_rxns=find_ex_rxns(model);
organic=find_organic_ex_rxns(model,ex_rxns);
model=changeRxnBounds(model,organic,-10,'l');
model.b=zeros(2774,1);
model=changeObjective(model,'biomassRxn');
optimizeCbModel(model)

confidenceScores=[];
 for i=1:length(model.confidenceScores)
     if isempty(model.confidenceScores{i})
         confidenceScores(i)=0;
     else
         confidenceScores(i)=str2num(sprintf('%.6f',model.confidenceScores{i}));
     end
 end
 

%section 1
evidence=tdfread('Intestine_bulk_expr_Section_1.txt','\t');
G=evidence.genes;
U=evidence.state;
G=cellstr(num2str(G));
G=strtrim(G);
C_H_genes=[];

%Add NAMPT gene association%
model=removeRxns(model,'NMNS');
model = addReaction(model,'NMMS',{'h[c]','ncam[c]','prpp[c]','ppi[c]','nmn[c]'},[-1 -1 -1 1 1],false,-1000,1000,0,'','59027.1');


[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, 0, C_H_genes, 2);
save('mouse_intestine_section1_model_mean_new_conn_changed_to_nansum_force_o2t.mat','PM','GM','C','NC','Z','model_C','pruneTime','cRes');


%section 5
load('MouseModel_New.mat')
model = MouseModel_New;
model=addReaction(model,'biomassRxn','0.60552 ala-L[c] + 0.030122 amp[c] + 0.39141 arg-L[c] + 0.33306 asn-L[c] + 0.26882 asp-L[c] + 30 atp[c] + 0.026219 chsterol[c] + 0.0019665 clpn_hs[c] + 0.030122 cmp[c] + 0.02 crn[c] + 0.16301 cys-L[c] + 0.003769 dag_hs[c] + 0.013565 damp[c] + 0.009074 dcmp[c] + 0.009074 dgmp[c] + 0.013565 dtmp[c] + 0.32832 gln-L[c] + 0.39909 glu-L[c] + 0.51919 gly[c] + 0.34069 glygn2[c] + 0.050294 gmp[c] + 0.01 gthrd[c] + 30 h2o[c] + 0.0074453 hdca[c] + 0.0029781 hdcea[c] + 0.14219 his-L[c] + 0.32929 ile-L[c] + 0.56489 leu-L[c] + 0.0032774 lpchol_hs[c] + 0.57111 lys-L[c] + 0.003769 mag_hs[c] + 0.13315 met-L[c] + 0.003 nad[c] + 0.0003 nadp[c] + 0.0044672 ocdca[c] + 0.025224 ocdcea[c] + 0.030152 pa_hs[c] + 0.0098323 pail_hs[c] + 0.0025 paps[c] + 0.047195 pchol_hs[c] + 0.044245 pe_hs[c] + 0.17578 phe-L[c] + 0.24458 pro-L[c] + 0.015076 ps_hs[c] + 0.0014 ptrc[c] + 4.8648e-05 q10[m] + 0.3973 ser-L[c] + 0.0091768 sphmyln_hs[c] + 0.016 spmd[c] + 0.006 sprm[c] + 7e-05 thbpt[c] + 0.38916 thr-L[c] + 0.041207 trp-L[c] + 0.15451 tyr-L[c] + 0.056957 ump[c] + 0.39711 val-L[c] + 0.022286 xolest_hs[c]  -> 30 adp[c] + 30 h[c] + 30 pi[c]');
model=addReaction(model,'IPDPxc','ipdp[x] 	->	ipdp[c]');
model=addReaction(model,'TYRtm','tyr-L[c] 	->	tyr-L[m]');
model=addReaction(model,'Tyr-ggnt','Tyr-ggn[e] 	->	Tyr-ggn[c]');
model=addReaction(model,'lystopeplys','lys-L[e] 	->	peplys[e]');
model=addReaction(model,'pepslys_sink','pepslys[r] 	->');
model = update_mouse_model(model);
ex_rxns=find_ex_rxns(model);
organic=find_organic_ex_rxns(model,ex_rxns);
model=changeRxnBounds(model,organic,-10,'l');
%model.b=zeros(2766,1);
%model=changeObjective(model,'biomassRxn');
optimizeCbModel(model)

confidenceScores=[];
 for i=1:length(model.confidenceScores)
     if isempty(model.confidenceScores{i})
         confidenceScores(i)=0;
     else
         confidenceScores(i)=str2num(sprintf('%.6f',model.confidenceScores{i}));
     end
 end
 
 

evidence=tdfread('Intestine_bulk_expr_Section_5.txt','\t');
G=evidence.genes;
U=evidence.state;
G=cellstr(num2str(G));
G=strtrim(G);
C_H_genes=[];

%Add NAMPT gene association%
model=removeRxns(model,'NMNS');
model = addReaction(model,'NMMS',{'h[c]','ncam[c]','prpp[c]','ppi[c]','nmn[c]'},[-1 -1 -1 1 1],false,-1000,1000,0,'','10135.1');


[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, 0, C_H_genes, 2);
save('mouse_intestine_section5_model_mean_new_conn_changed_to_nansum_force_o2t.mat','PM','GM','C','NC','Z','model_C','pruneTime','cRes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%without adding metabolites
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initCobraToolbox
%load('mouseModel.mat')
load('MouseModel_New.mat')
model = MouseModel_New;
model = update_mouse_model(model);
%load('precursorMets.mat')
%save('precursorMets_mmu.mat', 'precursorMets')
ex_rxns=find_ex_rxns(model);
organic=find_organic_ex_rxns(model,ex_rxns);
model=changeRxnBounds(model,organic,-10,'l');
%model.b=zeros(2766,1);
%model=changeObjective(model,'biomassRxn');
optimizeCbModel(model)

confidenceScores=[];
 for i=1:length(model.confidenceScores)
     if isempty(model.confidenceScores{i})
         confidenceScores(i)=0;
     else
         confidenceScores(i)=str2num(sprintf('%.6f',model.confidenceScores{i}));
     end
 end
 

%section 1
evidence=tdfread('Intestine_bulk_expr_Section_1.txt','\t');
G=evidence.genes;
U=evidence.state;
G=cellstr(num2str(G));
G=strtrim(G);
C_H_genes=[];


[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, 0, C_H_genes, 2);
save('mouse_intestine_section1_model_mean_new_conn_changed_to_nansum_force_o2t_no_add.mat','PM','GM','C','NC','Z','model_C','pruneTime','cRes');


%section 5
load('MouseModel_New.mat')
model = MouseModel_New;
model = update_mouse_model(model);
ex_rxns=find_ex_rxns(model);
organic=find_organic_ex_rxns(model,ex_rxns);
model=changeRxnBounds(model,organic,-10,'l');
%model.b=zeros(2766,1);
%model=changeObjective(model,'biomassRxn');
optimizeCbModel(model)

confidenceScores=[];
 for i=1:length(model.confidenceScores)
     if isempty(model.confidenceScores{i})
         confidenceScores(i)=0;
     else
         confidenceScores(i)=str2num(sprintf('%.6f',model.confidenceScores{i}));
     end
 end
 
 
 
evidence=tdfread('Intestine_bulk_expr_Section_5.txt','\t');
G=evidence.genes;
U=evidence.state;
G=cellstr(num2str(G));
G=strtrim(G);
C_H_genes=[];

%Add NAMPT gene association%
model=removeRxns(model,'NMNS');
model = addReaction(model,'NMMS',{'h[c]','ncam[c]','prpp[c]','ppi[c]','nmn[c]'},[-1 -1 -1 1 1],false,-1000,1000,0,'','10135.1');


[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, 0, C_H_genes, 2);
save('mouse_intestine_section5_model_mean_new_conn_changed_to_nansum_force_o2t_no_add.mat','PM','GM','C','NC','Z','model_C','pruneTime','cRes');
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%after conversion to human gene ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initCobraToolbox
load('recon1.mat')
model=addReaction(model,'biomassRxn','0.60552 ala-L[c] + 0.030122 amp[c] + 0.39141 arg-L[c] + 0.33306 asn-L[c] + 0.26882 asp-L[c] + 30 atp[c] + 0.026219 chsterol[c] + 0.0019665 clpn_hs[c] + 0.030122 cmp[c] + 0.02 crn[c] + 0.16301 cys-L[c] + 0.003769 dag_hs[c] + 0.013565 damp[c] + 0.009074 dcmp[c] + 0.009074 dgmp[c] + 0.013565 dtmp[c] + 0.32832 gln-L[c] + 0.39909 glu-L[c] + 0.51919 gly[c] + 0.34069 glygn2[c] + 0.050294 gmp[c] + 0.01 gthrd[c] + 30 h2o[c] + 0.0074453 hdca[c] + 0.0029781 hdcea[c] + 0.14219 his-L[c] + 0.32929 ile-L[c] + 0.56489 leu-L[c] + 0.0032774 lpchol_hs[c] + 0.57111 lys-L[c] + 0.003769 mag_hs[c] + 0.13315 met-L[c] + 0.003 nad[c] + 0.0003 nadp[c] + 0.0044672 ocdca[c] + 0.025224 ocdcea[c] + 0.030152 pa_hs[c] + 0.0098323 pail_hs[c] + 0.0025 paps[c] + 0.047195 pchol_hs[c] + 0.044245 pe_hs[c] + 0.17578 phe-L[c] + 0.24458 pro-L[c] + 0.015076 ps_hs[c] + 0.0014 ptrc[c] + 4.8648e-05 q10[m] + 0.3973 ser-L[c] + 0.0091768 sphmyln_hs[c] + 0.016 spmd[c] + 0.006 sprm[c] + 7e-05 thbpt[c] + 0.38916 thr-L[c] + 0.041207 trp-L[c] + 0.15451 tyr-L[c] + 0.056957 ump[c] + 0.39711 val-L[c] + 0.022286 xolest_hs[c]  -> 30 adp[c] + 30 h[c] + 30 pi[c]');
model=addReaction(model,'IPDPxc','ipdp[x] 	->	ipdp[c]');
model=addReaction(model,'TYRtm','tyr-L[c] 	->	tyr-L[m]');
model=addReaction(model,'Tyr-ggnt','Tyr-ggn[e] 	->	Tyr-ggn[c]');
model=addReaction(model,'lystopeplys','lys-L[e] 	->	peplys[e]');
model=addReaction(model,'pepslys_sink','pepslys[r] 	->');
ex_rxns=find_ex_rxns(model);
organic=find_organic_ex_rxns(model,ex_rxns);
model=changeRxnBounds(model,organic,-10,'l');
model.b=zeros(2766,1);
model=changeObjective(model,'biomassRxn');
optimizeCbModel(model)

confidenceScores=[];
 for i=1:length(model.confidenceScores)
     if isempty(model.confidenceScores{i})
         confidenceScores(i)=0;
     else
         confidenceScores(i)=str2num(model.confidenceScores{i});
     end
 end
 
evidence=tdfread('Intestine_bulk_expr_mouse2humanSection_1.txt','\t');
G=evidence.genes;
U=evidence.state;
G=cellstr(num2str(G));
G=strtrim(G);
C_H_genes=[];

%Add NAMPT gene association%
model=removeRxns(model,'NMNS');
model = addReaction(model,'NMMS',{'h[c]','ncam[c]','prpp[c]','ppi[c]','nmn[c]'},[-1 -1 -1 1 1],false,-1000,1000,0,'','10135.1');


[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, 0, C_H_genes, 2);

save('mouse_intestine_section1_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human','PM','GM','C','NC','Z','model_C','pruneTime','cRes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%after conversion to human gene ID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initCobraToolbox
load('recon1.mat')
model=addReaction(model,'biomassRxn','0.60552 ala-L[c] + 0.030122 amp[c] + 0.39141 arg-L[c] + 0.33306 asn-L[c] + 0.26882 asp-L[c] + 30 atp[c] + 0.026219 chsterol[c] + 0.0019665 clpn_hs[c] + 0.030122 cmp[c] + 0.02 crn[c] + 0.16301 cys-L[c] + 0.003769 dag_hs[c] + 0.013565 damp[c] + 0.009074 dcmp[c] + 0.009074 dgmp[c] + 0.013565 dtmp[c] + 0.32832 gln-L[c] + 0.39909 glu-L[c] + 0.51919 gly[c] + 0.34069 glygn2[c] + 0.050294 gmp[c] + 0.01 gthrd[c] + 30 h2o[c] + 0.0074453 hdca[c] + 0.0029781 hdcea[c] + 0.14219 his-L[c] + 0.32929 ile-L[c] + 0.56489 leu-L[c] + 0.0032774 lpchol_hs[c] + 0.57111 lys-L[c] + 0.003769 mag_hs[c] + 0.13315 met-L[c] + 0.003 nad[c] + 0.0003 nadp[c] + 0.0044672 ocdca[c] + 0.025224 ocdcea[c] + 0.030152 pa_hs[c] + 0.0098323 pail_hs[c] + 0.0025 paps[c] + 0.047195 pchol_hs[c] + 0.044245 pe_hs[c] + 0.17578 phe-L[c] + 0.24458 pro-L[c] + 0.015076 ps_hs[c] + 0.0014 ptrc[c] + 4.8648e-05 q10[m] + 0.3973 ser-L[c] + 0.0091768 sphmyln_hs[c] + 0.016 spmd[c] + 0.006 sprm[c] + 7e-05 thbpt[c] + 0.38916 thr-L[c] + 0.041207 trp-L[c] + 0.15451 tyr-L[c] + 0.056957 ump[c] + 0.39711 val-L[c] + 0.022286 xolest_hs[c]  -> 30 adp[c] + 30 h[c] + 30 pi[c]');
model=addReaction(model,'IPDPxc','ipdp[x] 	->	ipdp[c]');
model=addReaction(model,'TYRtm','tyr-L[c] 	->	tyr-L[m]');
model=addReaction(model,'Tyr-ggnt','Tyr-ggn[e] 	->	Tyr-ggn[c]');
model=addReaction(model,'lystopeplys','lys-L[e] 	->	peplys[e]');
model=addReaction(model,'pepslys_sink','pepslys[r] 	->');
ex_rxns=find_ex_rxns(model);
organic=find_organic_ex_rxns(model,ex_rxns);
model=changeRxnBounds(model,organic,-10,'l');
model.b=zeros(2766,1);
model=changeObjective(model,'biomassRxn');
optimizeCbModel(model)

confidenceScores=[];
 for i=1:length(model.confidenceScores)
     if isempty(model.confidenceScores{i})
         confidenceScores(i)=0;
     else
         confidenceScores(i)=str2num(model.confidenceScores{i});
     end
 end
 
evidence=tdfread('Intestine_bulk_expr_mouse2human_relative_section1.txt','\t');
G=evidence.genes;
U=evidence.state;
G=cellstr(num2str(G));
G=strtrim(G);
C_H_genes=[];

%Add NAMPT gene association%
model=removeRxns(model,'NMNS');
model = addReaction(model,'NMMS',{'h[c]','ncam[c]','prpp[c]','ppi[c]','nmn[c]'},[-1 -1 -1 1 1],false,-1000,1000,0,'','10135.1');


[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, 0, C_H_genes, 2);

save('mouse_intestine_relative_section1_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human','PM','GM','C','NC','Z','model_C','pruneTime','cRes');


load('recon1.mat')
model=addReaction(model,'biomassRxn','0.60552 ala-L[c] + 0.030122 amp[c] + 0.39141 arg-L[c] + 0.33306 asn-L[c] + 0.26882 asp-L[c] + 30 atp[c] + 0.026219 chsterol[c] + 0.0019665 clpn_hs[c] + 0.030122 cmp[c] + 0.02 crn[c] + 0.16301 cys-L[c] + 0.003769 dag_hs[c] + 0.013565 damp[c] + 0.009074 dcmp[c] + 0.009074 dgmp[c] + 0.013565 dtmp[c] + 0.32832 gln-L[c] + 0.39909 glu-L[c] + 0.51919 gly[c] + 0.34069 glygn2[c] + 0.050294 gmp[c] + 0.01 gthrd[c] + 30 h2o[c] + 0.0074453 hdca[c] + 0.0029781 hdcea[c] + 0.14219 his-L[c] + 0.32929 ile-L[c] + 0.56489 leu-L[c] + 0.0032774 lpchol_hs[c] + 0.57111 lys-L[c] + 0.003769 mag_hs[c] + 0.13315 met-L[c] + 0.003 nad[c] + 0.0003 nadp[c] + 0.0044672 ocdca[c] + 0.025224 ocdcea[c] + 0.030152 pa_hs[c] + 0.0098323 pail_hs[c] + 0.0025 paps[c] + 0.047195 pchol_hs[c] + 0.044245 pe_hs[c] + 0.17578 phe-L[c] + 0.24458 pro-L[c] + 0.015076 ps_hs[c] + 0.0014 ptrc[c] + 4.8648e-05 q10[m] + 0.3973 ser-L[c] + 0.0091768 sphmyln_hs[c] + 0.016 spmd[c] + 0.006 sprm[c] + 7e-05 thbpt[c] + 0.38916 thr-L[c] + 0.041207 trp-L[c] + 0.15451 tyr-L[c] + 0.056957 ump[c] + 0.39711 val-L[c] + 0.022286 xolest_hs[c]  -> 30 adp[c] + 30 h[c] + 30 pi[c]');
model=addReaction(model,'IPDPxc','ipdp[x] 	->	ipdp[c]');
model=addReaction(model,'TYRtm','tyr-L[c] 	->	tyr-L[m]');
model=addReaction(model,'Tyr-ggnt','Tyr-ggn[e] 	->	Tyr-ggn[c]');
model=addReaction(model,'lystopeplys','lys-L[e] 	->	peplys[e]');
model=addReaction(model,'pepslys_sink','pepslys[r] 	->');
ex_rxns=find_ex_rxns(model);
organic=find_organic_ex_rxns(model,ex_rxns);
model=changeRxnBounds(model,organic,-10,'l');
model.b=zeros(2766,1);
model=changeObjective(model,'biomassRxn');
optimizeCbModel(model)

confidenceScores=[];
 for i=1:length(model.confidenceScores)
     if isempty(model.confidenceScores{i})
         confidenceScores(i)=0;
     else
         confidenceScores(i)=str2num(model.confidenceScores{i});
     end
 end
 
evidence=tdfread('Intestine_bulk_expr_mouse2human_relative_section5.txt','\t');
G=evidence.genes;
U=evidence.state;
G=cellstr(num2str(G));
G=strtrim(G);
C_H_genes=[];

%Add NAMPT gene association%
model=removeRxns(model,'NMNS');
model = addReaction(model,'NMMS',{'h[c]','ncam[c]','prpp[c]','ppi[c]','nmn[c]'},[-1 -1 -1 1 1],false,-1000,1000,0,'','10135.1');


[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = ...
    mcadre(model, G, U, confidenceScores, 0, C_H_genes, 2);

save('mouse_intestine_relative_section5_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human','PM','GM','C','NC','Z','model_C','pruneTime','cRes');
