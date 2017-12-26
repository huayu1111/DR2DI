function candiRankList = rls_kron_for_diseases(y,ka,kb,DrugIDs,candiDrugs,DiseaseIndex,sigma)
	% Kronecker product Regularized Least Squares for association prediction.
	[va,la] = eig(ka);
	[vb,lb] = eig(kb);
	
	l = kron(diag(lb)',diag(la));
	l = l ./ (l + sigma);
	
	m1 = va' * y * vb;
	m2 = m1 .* l;
	y2 = va * m2 * vb';
	
	candiValueList=zeros(length(candiDrugs),1);
	for i=1:length(candiDrugs)
		logicIndex=strcmp(candiDrugs{i},DrugIDs);
		matchIndex = find(logicIndex,1);
		candiValueList(i) = y2(matchIndex,DiseaseIndex);
	end
	candiRankList(:,1)=candiDrugs;
	candiRankList(:,3)=num2cell(candiValueList);
end