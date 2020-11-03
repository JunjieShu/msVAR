

function yejiguiyin(prtf_test_r, portfolio_return_path, rownames, folder, k, p)


independent_variables = [ones(size(prtf_test_r,1),1), prtf_test_r.Variables];

[b,~,r,~,~] = regress( portfolio_return_path , independent_variables);   
load = b(2:end)';


contribution_rate = b(2:end)' .*  prtf_test_r.Variables   ./  portfolio_return_path;


mean_contribution_rate = mean(contribution_rate, 1);


h = figure;


theta = 0: 2*pi/(size(prtf_test_r,2)) : 2*pi - 2*pi/(size(prtf_test_r,2));
polarplot(theta, mean_contribution_rate, '-o')
rlim([min( mean_contribution_rate), max(mean_contribution_rate)])

thetaticks(rad2deg(theta))
thetaticklabels(rownames)

saveas(h,[folder,'/yejiguiyin(', num2str(k),',',num2str(p) ,').png'])
close(h)
