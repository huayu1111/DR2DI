function [ABCVCandiRankListDiseases] = ab_initio_cross_validation_for_diseases(datapath)
load([datapath,'dr2dinet.mat'])
load([datapath,'DrugIDs.mat']);
load([datapath,'DiseaseIDs.mat']);
load([datapath,'Drug2DiseaseRelations.mat']);
load([datapath,'DrugKernel.mat']);
load([datapath,'DiseaseKernel.mat']);
FinalCandiRankList=cell(0,0);
for i=1:length(DiseaseIDs) %#ok<USENS>
    logicIndex=strcmp(DiseaseIDs{i},Drug2DiseaseRelations(:,2)); %#ok<NODEF>
    DrugBin=Drug2DiseaseRelations(logicIndex,1);
    Tempdr2dinet = dr2dinet;
    for j=1:length(DrugBin)
        RowIndex=find(strcmp(DrugBin{j},DrugIDs),1);
        Tempdr2dinet(RowIndex,i)=0;
	end
	CandiRankList=rls_kron_for_diseases(Tempdr2dinet,drug_kernel_matrix,disease_kernel_matrix,DrugIDs,DrugIDs,i,1);
	CandiRankList(:,2)=DiseaseIDs(i);
    FinalCandiRankList=[FinalCandiRankList;CandiRankList]; %#ok<AGROW>
end
[~,IndexRow]=sort(cell2mat(FinalCandiRankList(:,3)),'descend');
FinalCandiRankList=FinalCandiRankList(IndexRow,:);
GoldData=strcat(Drug2DiseaseRelations(:,1),Drug2DiseaseRelations(:,2));
for i=1:length(FinalCandiRankList)
    BridgePair=strcat(FinalCandiRankList{i,1},FinalCandiRankList{i,2});
    if(~isempty(find(strcmp(BridgePair,GoldData),1)))
		% disp('Ok')
        FinalCandiRankList(i,4)={1};
        FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    else
        FinalCandiRankList(i,4)={0};
        FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    end
end
ABCVCandiRankListDiseases=FinalCandiRankList;
end