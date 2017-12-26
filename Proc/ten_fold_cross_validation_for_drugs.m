function [TenCVCandiRankListDrugs] = ten_fold_cross_validation_for_drugs(datapath)
load([datapath,'dr2dinet.mat'])
load([datapath,'DrugIDs.mat']);
load([datapath,'DiseaseIDs.mat']);
load([datapath,'Drug2DiseaseRelations.mat']);
load([datapath,'DrugKernel.mat']);
load([datapath,'DiseaseKernel.mat']);
FinalCandiRankList=cell(0,0);
for i=1:10 %#ok<USENS>
    if(i==10)
        Tempdr2dinet = dr2dinet;
        for j=(((i-1)*floor(length(DrugIDs)/10))+1):length(DrugIDs)
            logicIndex=strcmp(DrugIDs{j},Drug2DiseaseRelations(:,1)); %#ok<NODEF>
            DiseaseBin=Drug2DiseaseRelations(logicIndex,2);
            for k=1:length(DiseaseBin)
                ColIndex=find(strcmp(DiseaseBin{k},DiseaseIDs),1);
                Tempdr2dinet(j,ColIndex)=0;
            end
        end
        qDrugIndexes = (((i-1)*floor(length(DrugIDs)/10))+1):(i*floor(length(DrugIDs)/10));
        CandiRankList=rls_kron_for_ten_fold(Tempdr2dinet,drug_kernel_matrix,disease_kernel_matrix,DiseaseIDs,qDrugIndexes,DrugIDs,Drug2DiseaseRelations,1);
        FinalCandiRankList=[FinalCandiRankList;CandiRankList]; %#ok<AGROW>    
        
    else
        Tempdr2dinet = dr2dinet;
        for j=(((i-1)*floor(length(DrugIDs)/10))+1):(i*floor(length(DrugIDs)/10))
            logicIndex=strcmp(DrugIDs{j},Drug2DiseaseRelations(:,1)); %#ok<NODEF>
            DiseaseBin=Drug2DiseaseRelations(logicIndex,2);
            for k=1:length(DiseaseBin)
                ColIndex=find(strcmp(DiseaseBin{k},DiseaseIDs),1);
                Tempdr2dinet(j,ColIndex)=0;
            end
        end
        qDrugIndexes = (((i-1)*floor(length(DrugIDs)/10))+1):(i*floor(length(DrugIDs)/10));
        CandiRankList=rls_kron_for_ten_fold(Tempdr2dinet,drug_kernel_matrix,disease_kernel_matrix,DiseaseIDs,qDrugIndexes,DrugIDs,Drug2DiseaseRelations,1);
        FinalCandiRankList=[FinalCandiRankList;CandiRankList]; %#ok<AGROW>
    end
end
[~,IndexRow]=sort(cell2mat(FinalCandiRankList(:,3)),'descend');
FinalCandiRankList=FinalCandiRankList(IndexRow,:);
GoldData=strcat(Drug2DiseaseRelations(:,1),Drug2DiseaseRelations(:,2));
for i=1:length(FinalCandiRankList)
    BridgePair=strcat(FinalCandiRankList{i,1},FinalCandiRankList{i,2});
    if(~isempty(find(strcmp(BridgePair,GoldData),1)))
        FinalCandiRankList(i,4)={1};
        FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    else
        FinalCandiRankList(i,4)={0};
        FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    end
end
TenCVCandiRankListDrugs=FinalCandiRankList;
end