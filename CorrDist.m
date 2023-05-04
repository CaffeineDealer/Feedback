function [r p ch] = CorrDist(x)





for i = 1:size(x,2)
    x(i).df = x(i).df ./ max(x(i).df);
end

r = cell(size(x,2),1);
ch = r;
p = r;
for i = 1:size(x,2)
    for j = 1:size(x(i).df,2)
        if j ~= size(x(i).df,2)
            [R P] = corrcoef(x(i).df(:,j),x(i).df(:,j+1));
            ch{i} = [ch{i} x(i).ch(j+1) - x(i).ch(j)];
            r{i} = [r{i} R(1,2)];
            p{i} = [p{i} P(1,2)];
        end
    end
end


