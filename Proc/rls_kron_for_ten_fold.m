function finalCandiRankList = rls_kron_for_ten_fold(y,ka,kb,DiseaseIDs,DrugIndexes,DrugIDs,Drug2DiseaseRelations,sigma)
	% Kronecker product Regularized Least Squares for association prediction.
	
	[va,la] = eig(ka);
	[vb,lb] = eig(kb);
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	m1 = va' * y * vb;
	m2 = m1 .* l;
	y2 = va * m2 * vb';
	
	finalCandiRankList = cell(0,0);
	for i=1:length(DrugIndexes)
        logicIndexes=strcmp(DrugIDs{DrugIndexes(i)},Drug2DiseaseRelations(:,1)); 
		DiseaseBin=Drug2DiseaseRelations(logicIndexes,2);
        NonDiseaseBin = setdiff(DiseaseIDs,DiseaseBin);
        randomNonIndexes = randperm(length(NonDiseaseBin),length(DiseaseBin));
        candiDiseases=[DiseaseBin;NonDiseaseBin(randomNonIndexes)];
        candiRankList=cell(0,0);
		valueList=zeros(length(candiDiseases),1);
		for j=1:length(candiDiseases)
			logicIndex=strcmp(candiDiseases{j},DiseaseIDs);
			matchIndex = find(logicIndex,1);
			valueList(j) = y2(DrugIndexes(i),matchIndex);
		end
		candiRankList(:,2)=candiDiseases;
		candiRankList(:,3)=num2cell(valueList);
		candiRankList(:,1)=DrugIDs(DrugIndexes(i));
		finalCandiRankList=[finalCandiRankList;candiRankList];
	end
end
