function figures_states_on_assets(col, k,p, timetable_Returns, step,folder,markersize)
dates = timetable_Returns.Time;

new_folder_h = erase(join([folder, '/state_on_asset_return'])," ");
if ~exist(new_folder_h, 'dir')
       mkdir(new_folder_h)
end

new_folder_g = erase(join([folder, '/state_on_asset_CSI'])," ");
if ~exist(new_folder_g, 'dir')
       mkdir(new_folder_g)
end



n=size(timetable_Returns,2);
x0=10;
y0=10;
width=1500;
height=300*(n*2);

Markers = {'>','<','o','v','d','^','s','p','h','+','*','x'};
Colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30',...
               '#4DBEEE', '#A2142F' , 'r', 'g', 'm'};

%% state on asset return
h=figure('visible','off');
set(h,'position',[x0,y0,width,height])

for i = 1 : n
    subplot(n,1,i)
    legend_ = string([]); 
    hold on
    
    % plot asset
    plot(dates ,  timetable_Returns{: ,i}, ...
                    'LineStyle','-','LineWidth',1, 'Color','black')
    legend_=[legend_, {[timetable_Returns.Properties.VariableNames(i)]}];
    
    % plot states
    for j = 1 : k
        y = nan(size(col,1),1);
        temp = timetable_Returns{: ,i};
        y(col==j) = temp(col==j);
        
        plot(dates ,y, Markers{j},...
                    'MarkerSize',markersize,...
                    'MarkerEdgeColor',Colors{j},...
                    'MarkerFaceColor',Colors{j})          
        legend_ = [legend_, {['state', int2str(j),...
                    ' (', num2str(mean(y,'omitnan'),2), ', ', ...
                    num2str(var(y,'omitnan'),2),')']}]; 
                   
        plot(dates ,mean(y,'omitnan')*ones(size(col,1),1),...
                    '-.','Color',Colors{j},...
                    'LineWidth' ,1.5)    
                
        legend_ = [legend_, {['mean', int2str(j)]}]; 
    end
    title([timetable_Returns.Properties.VariableNames(i)])
    lgd = legend(legend_,'location','best');
    legend('boxoff');
    lgd.NumColumns = k+1;
    
    hold off
    grid on
   
end

sgtitle('state on asset return')


filename = erase(join([new_folder_h, '/(',  num2str(k),',',num2str(p)  ,').png'])," ");
saveas(h, filename)
close(h)




%% state on asset CSI
g=figure('visible','off');
set(g,'position',[x0,y0,width,height])

train_data_CSI = array2timetable(cumprod(exp(timetable_Returns.Variables/(250/step))), ...
                                   'RowTimes', timetable_Returns.Properties.RowTimes,...
                                   'VariableNames', timetable_Returns.Properties.VariableNames);

for i = 1 : n
    subplot(n,1,i)
    legend_ = string([]); 
    hold on
    plot(dates ,  train_data_CSI{: ,i}, ...
                    'LineStyle','-','LineWidth',2, 'Color','black')
    legend_=[legend_, {[train_data_CSI.Properties.VariableNames(i)]}];
    for j = 1 : k
        y = nan(size(col,1),1);
        temp = train_data_CSI{: ,i};
        y(col==j) = temp(col==j);
        
        plot(dates ,y, Markers{j},...
                    'MarkerSize',markersize,...
                    'MarkerEdgeColor',Colors{j},...
                    'MarkerFaceColor',Colors{j})
                
        legend_ = [legend_, {['state', int2str(j)]}];
    end
    title([train_data_CSI.Properties.VariableNames(i)])
    
    lgd = legend(legend_,'location','best');
    legend('boxoff');
    
    hold off
end


sgtitle('state on asset CSI')



%%
filename = erase(join([new_folder_g, '/(',  num2str(k),',',num2str(p)  ,').png'])," ");
saveas(g, filename)
close(g)


                            
                            