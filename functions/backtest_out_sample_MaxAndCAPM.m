
% Mdl = EstMdl

function [portfolio_weights, portfolio_return_path, msVAR_portfolio_annulized_return, ...
           equal_portfolio_annulized_return, excess_return_than_each_assets,...
           forecasted_SS]  ...
    = backtest_out_sample_MaxAndCAPM(Mdl, train_r, test_r,  prtf_train_r, prtf_test_r,prtf_train_rf ,prtf_test_rf,rownames)
% if no volatility control, set ctrl_vol = Inf

% reference:
% 
% msVAR/forecast function

% 
% AssetList = prtf_train_r.Properties.VariableNames;


% create folder
AssetList =  rownames;




P = Mdl.Switch.P;
k = Mdl.NumStates;
p = max(Mdl.P);
n = Mdl.NumSeries;

T = size(prtf_test_r, 1); % outsample length

r = [train_r.Variables; test_r.Variables];
prtf_r = [prtf_train_r.Variables; prtf_test_r.Variables];
prtf_rf = [prtf_train_rf.Variables; prtf_test_rf.Variables];


portfolio_weights = nan(size(prtf_test_r,1) , size(prtf_test_r,2)+1       );
portfolio_return = nan(size(prtf_test_r, 1), 1);

forecasted_SS = nan(size(prtf_test_r, 1), k);


%% 
pbar = CmdLineProgressBar('Out_sample Backtest...');
for t = 1 : T
    pbar.print(t, T)
    
    tt = t + size(prtf_train_r, 1); 
    
    % Y is used to find out the implied states
    Y = r(t:tt-1,:);
    prtf_Y = prtf_r(t:tt-1,:);
    

    [EstMdl,SS,~] = MSIAH(k, p, Y);
    forecasted_states = implied_states(SS);
    
    
    % Initial state probabilities
    S0 = SS(end,:);
    S0 = S0/sum(S0); % Normalize distribution
    S = S0(:)';
    
    % Begin forecast (Compute optimal forecast)
    S = S * EstMdl.Switch.P; % state probability of next period
    forecasted_SS(t,:) = S;
    
    
    
    weight = zeros(size(prtf_test_r,2)+1,  size(Mdl.Switch.P,2));
    
    for k = 1 : size(Mdl.Switch.P,2)

        AssetMean = mean(prtf_Y(forecasted_states==k, :), 1);
        AssetCovar = cov(prtf_Y(forecasted_states==k, :), 1);
        
        % if there is no such state, we put all the money into the test_rf
        if sum(isnan(AssetMean)) > 0
            weight(end,k) = 1;
            continue
        end

        % Create a Portfolio Object
        prft = Portfolio('AssetList',AssetList,'RiskFreeRate',test_rf{t,:});
        prft = setAssetMoments(prft,AssetMean,AssetCovar);
        prft = setDefaultConstraints(prft); % fully-invested, long-only portfolios

        % Maximize the Sharpe Ratio
        try   % sometimes the algorithm of matlab may have some error
            swgt = estimateMaxSharpeRatio(p);
        catch
            weight(end,k) = 1;
            continue
        end
        
        [srsk,sret] = estimatePortMoments(p,swgt);
        weight(1:end-1,k) = swgt;

   

    end
    
    % portfolio based on sharpe ratio and SS
    portfolio_weights(t,:) = S * weight';
    portfolio_return(t) = [ori_Y(tt,:), test_rf{t,:}] * weight * S';
  

end


portfolio_weights = array2timetable(portfolio_weights,...
                    'RowTimes', prtf_test_r.Time,...
                    'VariableNames', [prtf_test_r.Properties.VariableNames, {'Deposit'}] );




    equal_weight = ones(size(weight,1),1) / size(weight,1);
    equal_portfolio_return = [prtf_test_r, test_rf{:,:}] * equal_weight;
    best_portfolio_return = max([prtf_test_r.Variables,test_rf{:,:}],[],2);




portfolio_return_path = portfolio_return;
msVAR_portfolio_annulized_return = exp(sum(portfolio_return/12)) ^ (1/t*12) - 1;
equal_portfolio_annulized_return = exp(sum(equal_portfolio_return/12)) ^ (1/t*12) - 1;
% best_portfolio_annulized_return = exp(sum(best_portfolio_return/12)) ^ (1/t*12) - 1;

excess_return_than_each_assets = msVAR_portfolio_annulized_return - (exp(sum(prtf_test_r.Variables/12)) .^ (1/t*12) - 1);




% figure
% bar(prtf_test_r.Time, portfolio_weights.Variables,'stacked')
% legend(portfolio_weights.Properties.VariableNames)



%{

    pwgt = estimateFrontier(p,20);
    [prsk,pret] = estimatePortMoments(p,pwgt);
    q = setBudget(p, 0, 1);
    qwgt = estimateFrontier(q,20);
    [qrsk,qret] = estimatePortMoments(q,qwgt);            

    portfolioexamples_plot('E', ...
        {'line', prsk, pret}, ...
        {'line', qrsk, qret, [], [], 1}, ...
        {'scatter', srsk, sret, {'Sharpe'}}, ...
        {'scatter', vrsk, vret, {'ctrl'}}, ...
        {'scatter', sqrt(diag(p.AssetCovar)), p.AssetMean, p.AssetList, '.r'});


%}










