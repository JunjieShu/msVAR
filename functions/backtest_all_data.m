% Mdl = EstMdl


function  [real_SS, real_P] = backtest_all_data(Mdl, prtf_train_r, prtf_test_r)

k = Mdl.NumStates;
p = max(Mdl.P);


% calculate the SS of testing period using all data including train and test period


[real_Mdl,real_SS] = MSIAH(k, p, [prtf_train_r.Variables; prtf_test_r.Variables]);

real_P = real_Mdl.Switch.P;

real_SS = array2timetable(real_SS, 'RowTimes', [prtf_train_r.Time; prtf_test_r.Time],...
                        'VariableNames','k='+Mdl.StateNames);



end