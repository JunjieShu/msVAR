function figure_ss(k,p,SS, dates, folder)





new_folder = erase(join([folder, '/smooth_prob_and_stats'])," ");
if ~exist(new_folder, 'dir')
       mkdir(new_folder)
end


h=figure('visible','off');
x0=10;
y0=10;
width=2000;
height=300*(k+1);
set(h,'position',[x0,y0,width,height])


%% Smooth prob.
for i = 1 : k
    subplot(k+1,1,i)
    plot(dates, SS(:,i),'-','LineWidth',2)
    title(['state', int2str(i),' prob.' ])
end


%% impled states
subplot(k+1,1,i+1)
[impled_states,~]= ind2sub(size(SS'), find(SS'==max(SS',[],1)));
plot(dates, impled_states,'r','LineWidth',2)
yticks(1:k)
title('SS impied states')
grid on

%%
filename = erase(join([new_folder, '/ss_prob(',  num2str(k),',',num2str(p)  ,').png'])," ");
saveas(h, filename)




                            
                            