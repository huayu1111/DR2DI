function sim_mat = cal_sim()
	substructure = importdata('drug_substructure_mat.txt');
    target = importdata('drug_target_mat.txt');
    sider = importdata('drug_sider_mat.txt');
    S = size(substructure.data);
    substructure_freq = sum(substructure.data,1);
    substructure_std = std(substructure_freq);
    substructure_wei = exp(-(substructure_freq.^2)./((10*substructure_std)^2));
    target_freq = sum(target.data,1);
    target_std = std(target_freq);
    target_wei = exp(-(target_freq.^2)./((10*target_std)^2));
    sider_freq = sum(sider.data,1);
    sider_std = std(sider_freq);
    sider_wei = exp(-(sider_freq.^2)./((10*sider_std)^2));
    for k = 1:S(1)
        for j = (k+1):S(1)
            substructure_numerator = sum((substructure.data(k,:)).*(substructure.data(j,:)).*substructure_wei);
            substructure_denominator = sqrt(sum(substructure.data(k,:).*substructure_wei))*sqrt(sum(substructure.data(j,:).*substructure_wei));
            substructure_value = substructure_numerator/substructure_denominator;
            target_numerator = sum((target.data(k,:)).*(target.data(j,:)).*target_wei);
            target_denominator = sqrt(sum(target.data(k,:).*target_wei))*sqrt(sum(target.data(j,:).*target_wei));
            target_value = target_numerator/target_denominator;
            sider_numerator = sum((sider.data(k,:)).*(sider.data(j,:)).*sider_wei);
            sider_denominator = sqrt(sum(sider.data(k,:).*sider_wei))*sqrt(sum(sider.data(j,:).*sider_wei));
            sider_value = sider_numerator/sider_denominator;
            sim_mat(k,j) = max([substructure_value target_value sider_value]);
            sim_mat(j,k) = sim_mat(k,j);
        end
    end
end