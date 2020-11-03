function [abs_R] = abs_return(step, data_filepath, relative_return_filename, rownames     )


%% reading Primary Industry Index
filename = erase(join([data_filepath, relative_return_filename ,'.xlsx'])," ");


T = readtable(filename , 'Range', 'A10',...
                    'TreatAsEmpty',{'nan'}  );  
T=T(~any(ismissing(T),2), :);         

varTypes = repmat({'double'},1,size(T,2)-1 );
sz = [size(T,1), size(T,2)-1 ]; 

TT = timetable('Size', sz,'VariableTypes', varTypes, 'RowTimes', datetime(T.Var1) );
TT.Variables = table2array(T(:,2:end));

TT = sortrows(TT);
TT = rmmissing(TT,1);
TT.Properties.VariableNames = rownames;


temp = [];
datestrings = {};


for t = step+1 : step : size(TT,1)
    temp = [temp ; (250/step) * (log(TT{t,:}) - log(TT{t-step,:}))]; % annulized cc return
    datestrings = [datestrings , TT(t,:).Time];
    
end 

abs_R = array2timetable(temp,'RowTimes',datestrings, 'VariableNames',TT.Properties.VariableNames);













