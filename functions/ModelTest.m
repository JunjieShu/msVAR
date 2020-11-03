function  [JB, LR, LR_p_value] = ModelTest(EstMdl , SS, train_r,  z_max_p, z_max_l,  alpha, mctol)


%% Statistics
% disp('cacluating the z-stat...')

% for this formula, see Timmermann(2005) eqn 3.
Y = train_r.Variables;
VariableNames = train_r.Properties.VariableNames;

T = size(Y,1);
n = size(Y,2);
z = nan(T,n);
p = max(EstMdl.P);
k = EstMdl.NumStates;

% when we compute skew and kurt, the z_star cannot include inf data. To
% avoid this, we convert the 0 and 1 data in z to (1 / NumPaths)  and 
%  (NumPaths-1) / NumPaths

for n_ = 1 : n
    for t = p+1 : T
        
       z_t = 0;
        for i = 1 : k
            
            lag_term = 0;
            for j = 1 : p
                A_ij = cell2mat(EstMdl.Submodels(i).AR(j));
                lag_term = lag_term + A_ij(n_, n_) * Y(t-j ,n_);     
            end
            
            z_t = z_t + normcdf( (EstMdl.Submodels(i).Covariance(n_,n_))^(-1/2) *...
                (Y(t,n_) - EstMdl.Submodels(i).Constant(n_) -  lag_term) ) * SS(t, i);
            
        end
        
        z(t, n_) = z_t;

    end
end

z_star = norminv(z);


%% jbtest
% disp('executing the JBtest...')
% the null hypothesis that the data in vector x comes from a normal 
% distribution with an unknown mean and variance, using the Jarque-Bera 
% test. The alternative hypothesis is that it does not come from such a 
% distribution. 

% jbtest(x,alpha,mctol) returns a test decision based on a p-value computed
% using a Monte Carlo simulation with a maximum Monte Carlo standard error 
% less than or equal to mctol.

jb_h = nan(size(Y,2),1);
jb_p = nan(size(Y,2),1);
jb_stat = nan(size(Y,2),1);
Skew = nan(size(Y,2),1);
Kurt = nan(size(Y,2),1);
Average = nan(size(Y,2),1);
Sigma2 = nan(size(Y,2),1);

for i = 1 : size(Y,2)
    [jb_h(i),jb_p(i),jb_stat(i)] = jbtest(z_star(:,i), alpha, mctol);
    Skew(i) = skewness(z_star(:,i),0);
    Kurt(i) = kurtosis(z_star(:,i),0);
    Average(i) = mean(z_star(:,i),'omitnan');
    Sigma2(i) = var(z_star(:,i),0,'omitnan');
end


JB = [Average, Sigma2, Skew, Kurt, jb_stat, jb_p, jb_h]';


JB = array2table(JB,...
       'RowNames',{'ave','Var','Skew','Kurt', 'JB', 'p-value','h'},...
        'VariableNames',VariableNames);




%% LR test
% disp('executing the LR test...')
% see see Timmermann(2005) eqn 4 5 6
% null: good model ()
%

% Since for the OLS method, it is equivalent to maximize the likelihood
% function, so we firstly calculate the parameters by OLS, and then
% calculate equation5's maximized log-likelihood value. 
% reference: https://blog.csdn.net/solo_sky/article/details/47782657

% equation4
% OLS for equation5
% use i as loop indicater for z_max_l
% use j as indicator for z_max_p; this z_max_p is not the same as lag order p;
% use i as indcator for asset number: n
row_name = {};
p_ = [];
l_ = [];
LR = zeros(z_max_l*z_max_p, n);



for i = 1 : n          % asset class
    row_count = 0;             % reset row_count for next asset calculation
    
    for z_l = 1 : z_max_l         % lag order of test
        L =  - (T-z_l-p) /2*log(2*pi) -  sum(z_star( z_l+p+1:end, i).^2, 1)/2;  % Tim(2005)
        
        for z_p = 1 : z_max_p     % power          
            % construct independent variables
            independent_variables = NaN(T, z_l*z_p); 
            
                for t = z_l+p+1 : T     % sample size loop
                    for j = 1 : z_p
                        % for i = 1 : z_l % this for loop is replaced by
                        % matrix operation below
                        independent_variables(t,  (j-1)*z_l+1 : j*z_l) = (z_star(   t-z_l :t-1  ,i)').^j ;
                        % end
                    end
                end
            independent_variables = [ones(T,1),independent_variables];
            
            % regression
            [~,~,r,~,~] = regress(z_star(z_l+p+1:end, i), independent_variables(z_l+p+1:end, :));

            
            % maximized log-likelihood obtained from (5):
            row_count = row_count + 1;
            log_likelihood = (T-p-z_l)*log(1/sqrt(2*pi*var(r))) - 1/var(r)/2*sum(r.^2);
            LR(row_count, i) = -2*(L - log_likelihood);
            
            % save the row_name for table construction below the loops
            if i == 1
                row_name = [row_name, { ['p=',num2str(z_p),', l=',num2str(z_l) ]}];
                p_ = [p_, z_p];
                l_ = [l_, z_l];
            end

        end                      % end for z_p, the power in equation5
        
    end                          % end for lag order
    
end                              % end for asset class



% LR's p-value

LR_p_value = LR;

for i = 1 : n
    for row = 1 : size(LR,1)
        LR_p_value(row,i) = 1 - chi2cdf(LR(row, i),  p_(row)*l_(row) + 2 );
    end
end

LR_p_value = array2table( LR_p_value, ...
                'VariableNames',VariableNames,...
                'RowNames', row_name);
            
LR = array2table( LR, ...
                'VariableNames',VariableNames,...
                'RowNames', row_name); 


            
            
            
            


