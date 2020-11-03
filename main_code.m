tic

addpath('functions')

%% parameters

step = 5;

data_filepath = 'data/';

rf_filename = 'DepositInterestRate';
relative_return_filename = 'ZX-yjhyzs-daily'; 

rownames = string([ ...
        "石油石化", "煤炭", "有色金属", "电力及公用事业", "钢铁", "基础化工", ...
        "建筑", "建材", "轻工制造", "机械", "电力设备及新能源", "国防军工", ...
        "汽车", "商贸零售", "消费者服务", "家电", "纺织服装", "医药",...
        "食品饮料", "农林牧渔", "银行", "非银行金融", "房地产", ...
        "交通运输", "电子", "通信", "计算机", "传媒"]);




train_start_date = '2009-01-01';
train_end_date = '2018-12-31';
test_start_date = '2019-01-01';
test_end_date = '2019-12-31';

% save outputs in this folder
folder = ['results/step=',num2str(step), dbstack().name,'_' ,train_start_date,'~',train_end_date,';',test_start_date,'~',test_end_date];

if ~exist(folder, 'dir')
       mkdir(folder)
end



k_min = 3;
k_max = 6;
p_max = 2;
NumPaths = 1000;
alpha = 0.05;
mctol = 0.0001;
z_max_p = 2;
z_max_l = 2;
simu_num = 500; % numerical solution in function backtest_out_sample_Tim, same as in Timm 2005
gamma = 5;

% statistics are saved in this file
filename = erase(join([folder, '/', dbstack().name ,'_model_specification_test.xlsx'])," ");


%% convert to relative/absolute returns, considering the steps
% relative_R is used for pca and generate pca_data, and then train model;
% abs_R is used to invest when backtesting
[relative_R, market_R] = relative_return(step, data_filepath, relative_return_filename, rownames     );
abs_R = abs_return(step, data_filepath, relative_return_filename, rownames     );




%% dicrease the data dimension
X = table2array(relative_R);

% latent is the principal component variances, also  the eigenvalues of the covariance matrix of X
% score is the new coordinate after pca, where score = X * coeff
% explained is the Percentage of the total variance explained by each principal component,
[coeff,score,latent,tsquared,explained,mu] = pca(X);

clear result1
result1(1,:) = {'主元','贡献率', '累积贡献率'};
result1(2:size(X,2)+1, 1) = num2cell([1:size(X,2)]);
result1(2:size(X,2)+1, 2:3) = num2cell([explained, cumsum(explained)])



% 以元胞数组形式显示主成分表达式
clear result2
result2(:,1) = string(['industries', rownames]);

for i = 1 : 6
    result2(1,i+1) = {num2str(i)};
    result2(2:size(X,2)+1, i+1 ) = num2cell(coeff(:,i));
    
end

result2


%  plot scatter of three Principal components
figure
plot3(score(:,1),score(:,2), score(:,3), '+');
xlabel('1st Principal component');
ylabel('2nd Principal component');
zlabel('3rd Principal component');
grid on

clf




pca_data = array2timetable(score(:,1:5), 'RowTimes',relative_R.Time);



%% Data Process
% train model part
[train_r, test_r, train_rf, test_rf, train_dates, test_dates] = ...
    DataProcess(data_filepath, pca_data, rf_filename, ...
                     train_start_date, train_end_date,...
                    test_start_date, test_end_date );

                
% backtest part                
[prtf_train_r, prtf_test_r, prtf_train_rf, prtf_test_rf, ~, ~] = ...
    DataProcess(data_filepath, abs_R, rf_filename, ...
                     train_start_date, train_end_date,...
                    test_start_date, test_end_date );
                
                
% plot(table2array(pca_data))
                


% loop to find best parameters
for k = 4
   for p = 0 
    disp(['(   Below is ', 'k=',num2str(k),', p=', num2str(p),'   )']   )

%% estimate

% to avoid the error message when we set bad parameters (eg, to large k and
% p), which may cause numerical problem executing the program, we ues try
% catch statement to skip them.
try
    %%
    % there is a minor defference by the SS calculated by function MSIAH
    % (actually is function estimate) and function smooth. Since we use
    % smooth function in backtest, so here we also use smooth to
    % re-calculate SS. Of course we can abondon the first several SS to
    % make sure it's correct, but for now I skip this step. Perhaps I'll
    % refine it if needed in the future.
    [EstMdl,SS,logL,~,FS,S0] = MSIAH(k, p, train_r.Variables,'MaxIterations',999999);
    [SS,~,~,S0] = smooth(EstMdl, train_r.Variables);
    if p > 0
        SS = [repmat(S0,[p,1]);SS];
    end
    
    
    
    

    

catch
    continue
    
end

    %% test
    [JB, LR, LR_p_value] = ModelTest(EstMdl , SS, train_r,  z_max_p, z_max_l,  [], []);
    
    % if under such parameters, the model is hard to pass the test, then
    % just skip, no more figures
    
    if  sum(table2array(JB('p-value',:)) < 0.05) > (size(JB,2)/2)
        continue   
    end
    
    if sum(sum(table2array(LR_p_value) < 0.05)) > (size(LR_p_value,1) * size(LR_p_value,2) / 2)
        continue
    end
    
   
    
   
    %% Backtest in sample
    folder1 = [folder,'/1.in_sample'];
    
    if ~exist(folder1, 'dir')
       mkdir(folder1)
    end

    
    figure_ss(k,p,SS,train_dates,folder1);
    col = implied_states(SS);
    figures_states_on_assets(col, k,p, train_r, step,folder1, 5);


    %% Backtest out sample
    folder2 = [folder,'/2.out_sample'];
    
    if ~exist(folder2, 'dir')
       mkdir(folder2)
    end
    
%     [portfolio_weights, portfolio_return_path, msVAR_portfolio_annulized_return, ...
%            equal_portfolio_annulized_return, excess_return_than_each_assets,...
%            forecasted_SS]  ...
%     = backtest_out_sample_Max(EstMdl,train_r, test_r, prtf_train_r, prtf_test_r, prtf_test_rf,step);
% 

     [portfolio_weights, portfolio_return_path, msVAR_portfolio_annulized_return, ...
           equal_portfolio_annulized_return, excess_return_than_each_assets,...
           forecasted_SS]  ...
    = backtest_out_sample_Tim(EstMdl,train_r, test_r, train_r, test_r, prtf_train_rf,prtf_test_rf,step, simu_num, gamma);
    
    
    portfolio_weights_variables = portfolio_weights.Variables;
    pca_asset_weights = portfolio_weights_variables(:,1:end-1);
    
    
    weights = pca_asset_weights * coeff(:,1:5)';
    
    
    ori_weights = nan(size(prtf_test_r,1), size(prtf_test_r,2)+1 );
    
    ori_weights(:,end) = portfolio_weights.Deposit;
    
    ori_weights(:,1:end-1) =     weights ./ sum(weights,2) .* (1-portfolio_weights.Deposit);
    
    ori_weights = array2timetable(ori_weights,...
                    'RowTimes', prtf_test_r.Time,...
                    'VariableNames', [prtf_test_r.Properties.VariableNames, {'Deposit'}] );


  


    % culmulative return on portfolio in out-sample time
    h = figure;
    plot(test_dates,cumprod(exp(portfolio_return_path/(250/step))),'Linewidth',2)
    grid on
    saveas(h,[folder2,'/backtest_cul_ret_fig(', num2str(k),',',num2str(p) ,').png'])
    close(h)
    
    % portfolio weights in out-sample time
    h = figure;
    bar(table2array(ori_weights), 'stacked')
    grid on
    legend(ori_weights.Properties.VariableNames)
    saveas(h,[folder2,'/weights(', num2str(k),',',num2str(p) ,').png'])
    close(h)

    
    col_forecasted_SS = implied_states(forecasted_SS.Variables);
    

     figures_states_on_assets(col_forecasted_SS, k,p, test_r, step,folder2,10);

    disp(['msVAR:',num2str(100*msVAR_portfolio_annulized_return),...
    '% equal:',num2str(100*equal_portfolio_annulized_return),...
    ';;  % excess:',num2str(100*excess_return_than_each_assets),'  %'])




    %% portfolio performance
    marketReturn = market_R(test_dates, :).Variables;
    portReturn = portfolio_return_path;
    Cash = test_rf.Variables;
    
    [beta, alpha, sharpe_ratio, IR] = portfolio_performance(portReturn, marketReturn, Cash);
    
    
    %% Performance attribution
    yejiguiyin(prtf_test_r, portfolio_return_path, rownames, folder2, k, p)
    
   
  
    %% results output to excel
    sheet = ['MSIAH(', num2str(k), ',' , num2str(p),')' ] ;

    range_row = 1;
    writematrix('JB test', filename,  'Range', 'A1', 'Sheet', sheet   )
    range_row = range_row + 1;
    writetable(JB,filename, ...
        'WriteVariableNames',true, 'WriteRowNames',true ,...
        'Sheet', sheet,...
        'Range', ['A',num2str(range_row)]    )
    range_row = range_row + size(JB,1) + 2;

    writematrix('LR statistics', filename,  'Range', ['A',num2str(range_row)], 'Sheet', sheet   )
    range_row = range_row + 1;
    writetable(LR,filename, ...
        'WriteVariableNames',true, 'WriteRowNames',true ,...
        'Sheet', sheet,...
        'Range', ['A',num2str(range_row)])
    range_row = range_row + size(LR,1) + 2;

    writematrix('LR p-value', filename,  'Range', ['A',num2str(range_row)], 'Sheet', sheet   )
    range_row = range_row + 1;
    writetable(LR_p_value,filename, ...
        'WriteVariableNames',true, 'WriteRowNames',true ,...
        'Sheet', sheet,...
        'Range', ['A',num2str(range_row)])
    range_row = range_row + size(LR_p_value,1) + 2;
    
    range_row = range_row + 1;
    writematrix(excess_return_than_each_assets,filename, ...
        'Sheet', sheet,...
        'Range', ['A',num2str(range_row)])
   
    fprintf('\n')




   end
end


save([folder, '/all_variables(',  num2str(k),',',num2str(p)  , ').m'])

toc

