
% Mdl = EstMdl

function [portfolio_weights, portfolio_return_path, msVAR_portfolio_annulized_return, ...
       equal_portfolio_annulized_return, excess_return_than_each_assets,...
       forecasted_SS]  ...
    = backtest_out_sample_CAPM(Mdl, train_r, test_r, prtf_train_r, prtf_test_r,prtf_train_rf ,prtf_test_rf,step, simu_num, gamma)
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
prtf_rf = [prtf_train_rf.Variables; prtf_test_rf.Variables];



portfolio_weights = nan(size(prtf_test_r, 1), size(prtf_test_r, 2)+1);
portfolio_return = nan(size(prtf_test_r, 1), 1);
forecasted_SS = nan(size(prtf_test_r, 1), k);


%% 
pbar = CmdLineProgressBar('Out_sample Backtest...');
for t = 1 : T
    pbar.print(t, T)
    
    tt = t + size(prtf_train_r, 1); 
   
    
    % Y is used to find out the implied states
    Y = r(t:tt-1,:);
    [SS,~,~,S0] = smooth(Mdl,Y);
    if p > 0
        SS = [repmat(S0,[p,1]);SS];
    end  

   
    % initial SS to simulate data
    S0 = SS(end,:)'; % not same as S0 above
          
    
    % calculate portfolio weights on each assets
    %%%%%%% weight = zeros(size(prtf_test_r,2)+1,  size(Mdl.Switch.P,2));
    
  
    %% simulate
    % simulate to calculate utility to get portfolio weights and implied
    % states    
    % see Timm 2005 equ(15) or Barberis 2002 equ(6) 
     
    
    % R is the sum of risky asset return
    [temp_Y,~,statepath] = simulate(Mdl, T-t+1, 'S0', S0, 'NumPaths', simu_num );
    R = reshape( sum(temp_Y,1), simu_num,n  );
    
    % Rf is the sum of risky free asset return
    Rf = sum( ones(T-t+1,1)*prtf_rf(tt) );
    
    
    
    %%%% caluclate weights
    % w is used for function handle
    syms w
    for i=1:n
        syms (['w',num2str(i)]);
        eval(['w = [w;','w',num2str(i),'];']);
    end
    w(1) = [];
    
    % f is the Timm 2005 equ(15), aka expectation of the utility
    clear f  
    f = sum( ( (1-sum(w))*exp(Rf) + w.' * exp( Rf + R  )' ).^(1-gamma) / (1-gamma) ) / simu_num;
    f = matlabFunction(f);
    
    
    % numerical solution by using matlab function fmincon
    xx = 'x(1)';
    for i = 2 : n
        xx = [xx, ',x(',num2str(i),')'];
    end
    
    ff = ['@(x) -f(',xx,')'];
    
    % constrains
    w0 = ones(n,1)/n;
    A = ones(1,n);
    b = 1;
    ub = ones(n,1);
    lb = zeros(n,1);
    
    % sometimes the initial w0 is invalid for f1, then we use rand num to
    % avoid this problem
    while true
        try
            sol = fmincon( eval(ff), w0, A, b, [], [], lb, ub);
            break
        catch
            w0 = rand(n,1);
            w0 = w0 / sum(w0);
            continue
        end
    end
    
    
    portfolio_weights(t,:) = [sol' , 1 - sum(sol')];
    portfolio_return(t) = [prtf_r(tt,:), prtf_test_rf{t,:}] *  portfolio_weights(t,:)';
    
    %%%% calculate implied states SS
    
    for i = 1 : k
        forecasted_SS(t,i) =  (sum(statepath(1,:)==i)/simu_num) ;
    end
    
end


forecasted_SS = array2timetable(forecasted_SS, 'RowTimes', [prtf_test_r.Time],...
                        'VariableNames','k='+Mdl.StateNames);


portfolio_weights = array2timetable(portfolio_weights,...
                    'RowTimes', prtf_test_r.Time,...
                    'VariableNames', [prtf_test_r.Properties.VariableNames, {'Deposit'}] );



% calculate the equal weight for comparison
equal_weight = ones(size(portfolio_weights,2),1) / size(portfolio_weights,2);
equal_portfolio_return = [prtf_test_r.Variables, prtf_test_rf{:,:}] * equal_weight;
best_portfolio_return = max([prtf_test_r.Variables, prtf_test_rf{:,:}],[],2);




portfolio_return_path = portfolio_return;
msVAR_portfolio_annulized_return = exp(sum(portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;
equal_portfolio_annulized_return = exp(sum(equal_portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;
% best_portfolio_annulized_return = exp(sum(best_portfolio_return/(250/step))) ^ (1/t*(250/step)) - 1;


excess_return_than_each_assets = msVAR_portfolio_annulized_return - (exp(sum(prtf_test_r.Variables/(250/step))) .^ (1/t*(250/step)) - 1);











