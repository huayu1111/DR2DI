function [LOOCVCandiRankListDrugs] = leave_one_out_cross_validation_for_drugs(datapath)
load([datapath,'dr2dinet.mat'])
load([datapath,'DrugIDs.mat']);
load([datapath,'DiseaseIDs.mat']);
load([datapath,'Drug2DiseaseRelations.mat']);
load([datapath,'DrugKernel.mat']);
load([datapath,'DiseaseKernel.mat']);
FinalCandiRankList=cell(0,0);
for i=1:length(DrugIDs) %#ok<USENS>
    logicIndex=strcmp(DrugIDs{i},Drug2DiseaseRelations(:,1)); %#ok<NODEF>
    DiseaseBin=Drug2DiseaseRelations(logicIndex,2);
    TempCandiRankList=cell(0,0);
    Tempdr2dinet = dr2dinet;
    for j=1:length(DiseaseBin)
        ColIndex=find(strcmp(DiseaseBin{j},DiseaseIDs),1);
        Tempdr2dinet(i,ColIndex)=0;
        CandiDiseases=setdiff(DiseaseIDs,DiseaseBin);
        TempCandi=rls_kron_for_drugs(Tempdr2dinet,drug_kernel_matrix,disease_kernel_matrix,DiseaseIDs,[DiseaseBin{j};CandiDiseases],i,1);
        TempCandi(:,1)=DrugIDs(i);
        TempCandiRankList=[TempCandiRankList;TempCandi]; %#ok<AGROW>
    end
    DiseaseUniqs=unique(TempCandiRankList(:,2));
    CandiRankList=cell(0,0);
    for s=1:length(DiseaseUniqs)
        LogicIndex=strcmp(DiseaseUniqs{s},TempCandiRankList(:,2));
        ColIndex=find(LogicIndex);
        RankValue=mean(cell2mat(TempCandiRankList(LogicIndex,3)));
        CandiRankList=[CandiRankList;TempCandiRankList(ColIndex(1),:)]; %#ok<AGROW>
        CandiRankList(end,3)={RankValue};
    end
    FinalCandiRankList=[FinalCandiRankList;CandiRankList]; %#ok<AGROW>
end
[~,IndexRow]=sort(cell2mat(FinalCandiRankList(:,3)),'descend');
FinalCandiRankList=FinalCandiRankList(IndexRow,:);
GoldData=strcat(Drug2DiseaseRelations(:,1),Drug2DiseaseRelations(:,2));
for i=1:length(FinalCandiRankList)
    BridgePair=strcat(FinalCandiRankList{i,1},FinalCandiRankList{i,2});
    if(~isempty(find(strcmp(BridgePair,GoldData),1)));
		disp('Ok')
        FinalCandiRankList(i,4)={1};
        FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    else
        FinalCandiRankList(i,4)={0};
        FinalCandiRankList(i,5)={(length(FinalCandiRankList)-i)/length(FinalCandiRankList)};
    end
end
LOOCVCandiRankListDrugs=FinalCandiRankList;
end