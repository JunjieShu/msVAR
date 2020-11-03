
% Mdl = EstMdl;
% 


function [portfolio_return_path, msVAR_portfolio_annulized_return, ...
         equal_portfolio_annulized_return, best_portfolio_annulized_return]  ...
    = backtest_in_sample_max(Mdl, col, SS, train_Y, train_rf,step)
% if no volatility control, set ctrl_vol = Inf

% reference:
% 
% Portfolio Optimization Examples
% https://www.mathworks.com/help/finance/examples/portfolio-optimization-examples-1.html
%
% Portfolio Analysis
% https://www.mathworks.com/help/finance/portfolio-analysis.html
%
% Portfolio Optimization Against a Benchmark
% https://www.mathworks.com/help/finance/examples/portfolio-optimization-against-dow-benchmark.html



Data = train_Y.Variables;
T = size(train_Y, 1);

switch_P = Mdl.Switch.P;




% 
portfolio_return = nan(size(train_Y, 1), 1);

pbar = CmdLineProgressBar('In_sample Backtest...');
for t = 1 : T
    pbar.print(t,T)
    
    weight = zeros(size(train_Y,2)+1,  size(switch_P,2));
    for k = 1 : size(Mdl.Switch.P,2)
        AssetMean = mean(Data(col==k, :), 1);

        
        max_asset_mean = max(AssetMean);
        if max(AssetMean) < 0
            weight(end,k) = 1;
        else
            weight(find(AssetMean == max_asset_mean), k) = 1;
        end
    end
    
    
    portfolio_return(t) = [Data(t,:), train_rf{t,:}] * weight * SS(t,:)';
  
end



equal_weight = ones(size(weight,1),1) / size(weight,1);
equal_portfolio_return = [Data, train_rf{:,:}] * equal_weight;
best_portfolio_return = max([Data,train_rf{:,:}],[],2);




portfolio_return_path = portfolio_return;
msVAR_portfolio_annulized_return = exp(sum(portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;
equal_portfolio_annulized_return = exp(sum(equal_portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;
best_portfolio_annulized_return = exp(sum(best_portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;





%     pwgt = estimateFrontier(p,20);
%     [prsk,pret] = estimatePortMoments(p,pwgt);
%     q = setBudget(p, 0, 1);
%     qwgt = estimateFrontier(q,20);
%     [qrsk,qret] = estimatePortMoments(q,qwgt);            
% 
%     portfolioexamples_plot('E', ...
%         {'line', prsk, pret}, ...
%         {'line', qrsk, qret, [], [], 1}, ...
%         {'scatter', srsk, sret, {'Sharpe'}}, ...
%         {'scatter', vrsk, vret, {'ctrl'}}, ...
%         {'scatter', sqrt(diag(p.AssetCovar)), p.AssetMean, p.AssetList, '.r'});


















