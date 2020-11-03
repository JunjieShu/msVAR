
% Mdl = EstMdl

function [portfolio_weights, portfolio_return_path, msVAR_portfolio_annulized_return, ...
       equal_portfolio_annulized_return, excess_return_than_each_assets,...
       forecasted_SS]  ...
    = backtest_out_sample_Max(Mdl, train_r, test_r, prtf_train_r, prtf_test_r, prtf_test_rf,step)
% if no volatility control, set ctrl_vol = Inf

% reference:
% 
% msVAR/forecast function



P = Mdl.Switch.P;
k = Mdl.NumStates;
p = max(Mdl.P);
n = Mdl.NumSeries;

T = size(prtf_test_r, 1); % outsample length

r = [train_r.Variables; test_r.Variables];
prtf_r = [prtf_train_r.Variables; prtf_test_r.Variables];



portfolio_weights = nan(size(prtf_test_r, 1), size(prtf_test_r, 2)+1);
portfolio_return = nan(size(prtf_test_r, 1), 1);
forecasted_SS = nan(size(prtf_test_r, 1), k);


%% 
pbar = CmdLineProgressBar('Out_sample Backtest...');
for t = 1 : T
    pbar.print(t, T)
    
    tt = t + size(prtf_train_r, 1); 
    
    % original returns, used to calculate the final returns
    Data = prtf_r(t:tt-1,:);  
    

    % Y is used to find out the implied states
    Y = r(t:tt-1,:);
    [SS,~,~,S0] = smooth(Mdl,Y);
    if p > 0
        SS = [repmat(S0,[p,1]);SS];
    end  
    col = implied_states(SS);
    
    
    % calculate the state prob., by Tim 2005 (11)
    S0 = SS(end,:)'; % not same as S0 above
    
    eta = NaN(k,1);
    for k_ = 1 : k
        
        mu = Mdl.Submodels(k_).Constant ;
        for p_ = 1 : p
            mu = mu + cell2mat(Mdl.Submodels(k_).AR(p_)) * Y(end - p_, :)';
        end
              
        eta(k_) = mvnpdf( Y(end, : ), mu', Mdl.Submodels(k_).Covariance);
        
    end
    
    SS_next = ((S0'*P)' .* eta )   /   ( ((S0'*P)' .* eta )' * ones(k,1));
    forecasted_SS(t,:) = SS_next';
    
    
%     S0 = SS(end,:)'; % not same as S0 above 
%     SS_next = (S0'*P)' ;
%     forecasted_SS(t,:) = SS_next';
    
    
    % calculate portfolio weights on each assets
    weight = zeros(size(prtf_test_r,2)+1,  size(Mdl.Switch.P,2));
    
    for k_ = 1 : k
        
        AssetMean = mean(Data(col==k_, :), 1);
        
        % if there is no such state, we put all the money into the test_rf
        if sum(isnan(SS_next)) > 0
            weight(end,k_) = 1;
            continue
        end
        
        % if the return of the mean of all assets less than 0, we put all 
        % the money into the Deposit
        AssetMean_max = max(AssetMean);
        if max(AssetMean_max) < 0 || max(AssetMean_max) < prtf_test_rf{t,:}
            weight(end, k_) = 1;
        else
            weight(find(AssetMean == AssetMean_max), k_) = 1;
        end
   
    end
    
    
    portfolio_weights(t,:) = weight * SS_next;
    portfolio_return(t) = [prtf_r(tt,:), prtf_test_rf{t,:}] * weight * SS_next;
  

end


portfolio_weights = array2timetable(portfolio_weights,...
                    'RowTimes', prtf_test_r.Time,...
                    'VariableNames', [prtf_test_r.Properties.VariableNames, {'Deposit'}] );




equal_weight = ones(size(weight,1),1) / size(weight,1);
equal_portfolio_return = [prtf_test_r.Variables, prtf_test_rf{:,:}] * equal_weight;
best_portfolio_return = max([prtf_test_r.Variables, prtf_test_rf{:,:}],[],2);




portfolio_return_path = portfolio_return;
msVAR_portfolio_annulized_return = exp(sum(portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;
equal_portfolio_annulized_return = exp(sum(equal_portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;
% best_portfolio_annulized_return = exp(sum(best_portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;


excess_return_than_each_assets = msVAR_portfolio_annulized_return - (exp(sum(prtf_test_r.Variables/(250/step))) .^ (1/t*(250/step)) - 1);


forecasted_SS = array2timetable(forecasted_SS, 'RowTimes', [prtf_test_r.Time],...
                        'VariableNames','k='+Mdl.StateNames);











