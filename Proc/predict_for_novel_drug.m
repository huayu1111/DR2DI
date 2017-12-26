function [PredCandiRankList] = predict_for_novel_drug(datapath)
load([datapath,'dr2dinet.mat']);
load([datapath,'DrugIDs.mat']);
load([datapath,'DiseaseIDs.mat']);
load([datapath,'DrugKernel.mat']);
load([datapath,'DiseaseKernel.mat']);
PredCandiRankList=rls_kron_for_predict(dr2dinet,drug_kernel_matrix,disease_kernel_matrix,DiseaseIDs,DrugIDs,1);
end