function [train_r, test_r, train_rf, test_rf, train_dates, test_dates] = ...
    DataProcess(data_filepath, Return_data, rf_filename, ...
                     train_start_date, train_end_date,...
                    test_start_date, test_end_date )


% training and testing period
trainPrd = datetime([train_start_date; train_end_date],'InputFormat','yyyy-MM-dd');
testPrd = datetime([test_start_date; test_end_date],'InputFormat','yyyy-MM-dd');


% risk free rate
filename = erase(join([data_filepath, rf_filename,'.xlsx'])," ");
T = readtable(filename , 'Range', 'A10',...
                    'TreatAsEmpty',{'nan'});    
T=T(~any(ismissing(T),2), :);
rf = timetable(datetime(T.Var1) ,T.Var2);
rf = varfun(@(x) x/100, rf);

rf = sortrows(rf);
rf = rmmissing(rf,1);
rf.Properties.VariableNames = {'rf'};

% rf = retime(rf, 'daily','lastvalue');
% rf.Time.Format = 'uuuu-MM-D';





% Relative Returns 


train_r = Return_data(isbetween(Return_data.Time,trainPrd(1),trainPrd(2)),:);
test_r = Return_data(isbetween(Return_data.Time,testPrd(1),testPrd(2)),:);

train_rf = rf(train_r.Time,:);
test_rf = rf(test_r.Time,:);

train_dates = train_r.Time;
test_dates = test_r.Time;



