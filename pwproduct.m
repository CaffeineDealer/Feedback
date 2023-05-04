function [pwp] = pwproduct(x)

% Pairwise product function calculates multiplication between every pairs
% of a vector
% Written by Yavar Korkian on 2023/02/12
% Last Modified on 2023/02/12

y = allcomb(x,x);
% y(y(:,1)==y(:,2),:) = [];
y = y(:,1) .* y(:,2);
pwp = y;