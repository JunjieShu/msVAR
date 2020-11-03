function col = implied_states(ss)

[col,~]= ind2sub(size(ss'), find(ss'==max(ss',[],1)));