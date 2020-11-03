



function [beta, alpha, sharpe_ratio, IR] = portfolio_performance(portReturn, marketReturn, Cash)

beta = portbeta(portReturn, marketReturn);  %% my own function
alpha = portalpha(portReturn, marketReturn, Cash, 'capm'   );
sharpe_ratio = sharpe(portReturn, Cash);
[ IR , TE ] = inforatio ( portReturn , marketReturn );



% TRate = cumprod(exp(portfolio_return_path/(250/step)));
% [ MaxDD , MaxDDIndex ] = maxdrawdown ( TRate , 'arithmetic' )
% 
% clf
% figure
% 
% 
% plot(TRate)
% hold on
% plot( MaxDDIndex , TRate ( MaxDDIndex ), 'r-o' , 'MarkerSize' , 10)
% 
% 
% 
% 

