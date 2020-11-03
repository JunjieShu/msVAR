function [EstMdl,SS,logL,AIC,BIC,HQC] = MSIAH(k, p, Y, varargin)
%% Checking arguments
% In this scheme, we check the optional argument X (the predictor variable)
% Our model could include the predictor variables or not. In Timmermann
% (2005), it's the dividend yield. Since this argument is optional, we use
% the InputParser object to check its validation; only when b>0 and b is
% int this argument will be valid.

% Create an InputParser object

% defaultX = NaN;
DefaultIterationPlot = false;
DefaultMaxIterations = 1000;


parser = inputParser;
valid_Num = @(x) isnumeric(x);
valid_Logic = @(x) islogical(x);



addRequired(parser,'k',valid_Num);
addRequired(parser,'p',valid_Num);
addRequired(parser,'Y',valid_Num);

addParameter(parser,'IterationPlot',DefaultIterationPlot,valid_Logic);
addParameter(parser,'MaxIterations',DefaultMaxIterations,valid_Num);

parse(parser,k,p,Y,varargin{:});

k = parser.Results.k;
p = parser.Results.p;
Y = parser.Results.Y;

IterationPlot = parser.Results.IterationPlot;
MaxIterations = parser.Results.MaxIterations;

n = size(Y,2);


% no beta terms in models, all are endo variables
% b = 0;
para_num = k*(p*n^2 + n );
Obs_num = numel(Y);





%% Create (Partially Specified) Model for Estimation
P = nan(k);
mc = dtmc(P);
mdl = repmat(varm(n, p), k, 1);

% if b ~= 0
%     for i = 1 : k
%         mdl(i).Beta = NaN(n,b);
%     end
% end

Mdl = msVAR(mc, mdl);

%% Create Fully Specified Model Containing Initial Values
P0 = (1/k)*ones(k);
mc0 = dtmc(P0);
mdl0 = repmat(varm(n, p), k, 1);

Y_sort = sort(Y);
C = zeros(k,n);

% cvar = zeros(n,n,k);

% initiate AR
VAR_Mdl = varm(n,p);
Est_VAR_mdl = estimate(VAR_Mdl, Y);
AR0 = Est_VAR_mdl.AR;

for i = 1 : k
    mdl0(i).AR = AR0;
end

% % initiate Beta
% if b ~= 0
%     for i = 1 : k
%         mdl0(i).Beta = randn(n,b);
%     end
% end


% initiate Constant and Convariance

for i = 1 : k
    C(i,:) = mean(...
       Y_sort( 1+floor((i-1)*(1/k)*size(Y_sort,1)) : floor(i*(1/k)*size(Y_sort,1)), :) );
end


for i = 1 : k
    mdl0(i).Constant = C(i,:)';
    % mdl0(i).Covariance = cvar(:,:,i);
    mdl0(i).Covariance = eye(n);
    
   
end




Mdl0 = msVAR(mc0, mdl0);



%% Estimation
if IterationPlot
    figure
end



[EstMdl,SS,logL,~,~,S0] = estimate(Mdl,Mdl0,Y,'IterationPlot',IterationPlot,...
                            'MaxIterations', MaxIterations);


if p > 0
    % the first p Y are used by estimate function, so the first p time has
    % no SS, we use the initial S0 (steady-state distribution computed by 
    % DTMC/ASYMPTOTICS.)
    
    SS = [repmat(S0,[p,1]);SS];
end




AIC = -2*logL + 2*para_num;
BIC = -2*logL + para_num*log(Obs_num);
HQC = -2*logL + 2*para_num*log(log(Obs_num));


