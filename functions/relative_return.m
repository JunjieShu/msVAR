function [relative_R, market_R] = relative_return(step, data_filepath, relative_return_filename, rownames     )


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

R = array2timetable(temp,'RowTimes',datestrings, 'VariableNames',TT.Properties.VariableNames);


%% reading Market index 
filename = erase(join([data_filepath, 'HS300','.xlsx'])," ");
T = readtable(filename , 'Range', 'A10',...
                    'TreatAsEmpty',{'nan'}  );  
T=T(~any(ismissing(T),2), :);         

varTypes = repmat({'double'},1,size(T,2)-1 );
sz = [size(T,1), size(T,2)-1 ]; 

TTm = timetable('Size', sz,'VariableTypes', varTypes, 'RowTimes', datetime(T.Var1) );
TTm.Variables = table2array(T(:,2:end));

TTm = sortrows(TTm);
TTm = rmmissing(TTm,1);
TTm.Properties.VariableNames = "Market-index-szzhzs";

% make sure market index is in the range Primary Industry Index
TTm = TTm(   isbetween(TTm.Time, TT.Time(1), TT.Time(end)) , :  );


temp = [];
datestrings = {};


for t = step+1 : step : size(TTm,1)
    temp = [temp ; (250/step) * (log(TTm{t,:}) - log(TTm{t-step,:}))]; % annulized cc return
    datestrings = [datestrings , TTm(t,:).Time];
    
end 

market_R = array2timetable(temp,'RowTimes',datestrings, 'VariableNames',TTm.Properties.VariableNames);


%% relative return of Primary Industry Index to market index

if size(R,1) > size(market_R,1)
    R = R(market_R.Time, :);
elseif size(R,1) < size(market_R,1)
    market_R = market_R(R.Time, :);
end


relative_R = table2array(R) - table2array(market_R);

relative_R = array2timetable(relative_R, 'RowTimes',datestrings, 'VariableNames',rownames);









