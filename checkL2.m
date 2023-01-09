% function flag = checkL2(L2)
% 
% In:  L2   = QPack Retrieval Structure L2
%
% Out: flag = 1 = good, 0 = bad
%
%
% Dominik Scheiben, 2010-11-12
%
function flag = checkL2(L2)

p_ind = find(L2.species1_p/100>=0.01 & L2.species1_p/100<=10);
if isempty(find(L2.species1_x(p_ind)*1e6>12 | L2.species1_x(p_ind)*1e6<0, 1))
    flag = 1;
else
    flag = 0;
end