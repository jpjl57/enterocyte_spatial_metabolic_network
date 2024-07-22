initCobraToolbox
load('mouse_intestine_relative_section1_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat')
model1 = PM;

sec1 = findRxnsFromGenes(model1, '206358.1');
alt1 = findRxnsFromGenes(model1, '124935.1');
findRxnsFromGenes(model, '952.1');
findMetsFromRxns(model1, 'GUAPRT')
findRxnsFromMets(model1, 'nad[c]')
findGenesFromRxns(model, 'NMMS')

load('mouse_intestine_relative_section5_model_mean_new_conn_changed_to_nansum_force_o2t_mouse2human.mat');
model5 = PM;

sec5 = findRxnsFromGenes(model5, '2766.1');
alt5 = findRxnsFromGenes(model5, '124935.1');
findRxnsFromGenes(model, '952.1');
findMetsFromRxns(model5, 'SARDHm')
findRxnsFromMets(model5, 'ru5p-D[c]')
findGenesFromRxns(model, 'EX_ncam_e')
