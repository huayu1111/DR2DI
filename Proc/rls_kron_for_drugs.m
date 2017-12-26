function candiRankList = rls_kron_for_drugs(y,ka,kb,DiseaseIDs,candiDiseases,DrugIndex,sigma)
	% Kronecker product Regularized Least Squares for association prediction.
	
	[va,la] = eig(ka);
	[vb,lb] = eig(kb);
	
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	
	m1 = va' * y * vb;
	m2 = m1 .* l;
	y2 = va * m2 * vb';
	
	candiValueList=zeros(length(candiDiseases),1);
	for i=1:length(candiDiseases)
		logicIndex=strcmp(candiDiseases{i},DiseaseIDs);
		matchIndex = find(logicIndex,1);
		candiValueList(i) = y2(DrugIndex,matchIndex);
	end
	
	candiRankList(:,2)=candiDiseases;
	candiRankList(:,3)=num2cell(candiValueList);
end
