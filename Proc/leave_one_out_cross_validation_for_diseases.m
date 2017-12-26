function [LOOCVCandiRankListDiseases] = leave_one_out_cross_validation_for_diseases(datapath)
load([datapath,'dr2dinet.mat']);
load([datapath,'DrugIDs.mat']);
load([datapath,'DiseaseIDs.mat']);
load([datapath,'Drug2DiseaseRelations.mat']);
load([datapath,'DrugKernel.mat']);
load([datapath,'DiseaseKernel.mat']);
FinalCandiRankList=cell(0,0);
for i=1:length(DiseaseIDs) %#ok<USENS>
    logicIndex=strcmp(DiseaseIDs{i},Drug2DiseaseRelations(:,2)); %#ok<NODEF>
    DrugBin=Drug2DiseaseRelations(logicIndex,1);
    TempCandiRankList=cell(0,0);
    Tempdr2dinet = dr2dinet;
    for j=1:length(DrugBin)
        RowIndex=find(strcmp(DrugBin{j},DrugIDs),1);
        Tempdr2dinet(RowIndex,i)=0;
        CandiDrugs=setdiff(DrugIDs,DrugBin);
        TempCandi=rls_kron_for_diseases(Tempdr2dinet,drug_kernel_matrix,disease_kernel_matrix,DrugIDs,[DrugBin{j};CandiDrugs],i,1);
        TempCandi(:,2)=DiseaseIDs(i);
        TempCandiRankList=[TempCandiRankList;TempCandi]; %#ok<AGROW>
    end
    DrugUniqs=unique(TempCandiRankList(:,1));
    CandiRankList=cell(0,0);
    for s=1:length(DrugUniqs)
        LogicIndex=strcmp(DrugUniqs{s},TempCandiRankList(:,1));
        RowIndex=find(LogicIndex);
        RankValue=mean(cell2mat(TempCandiRankList(LogicIndex,3)));
        CandiRankList=[CandiRankList;TempCandiRankList(RowIndex(1),:)]; %#ok<AGROW>
        CandiRankList(end,3)={RankValue};
    end
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
LOOCVCandiRankListDiseases=FinalCandiRankList;
end