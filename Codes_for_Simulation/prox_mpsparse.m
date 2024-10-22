function [prox] = prox_mpsparse(x,m)

n = numel(x);
prox = x;
[~, sorted_ind] = sort(prox,'descend');
prox(sorted_ind(m+1:n)) = 0;
prox(prox<0) = 0;

end

