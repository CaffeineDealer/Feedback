function [w] = CorrDistMSTMT(x,y)





for i = 1:size(x,2)
    x(i).df = x(i).df ./ max(x(i).df);
end


for i = 1:size(x,2)
    for j = 1:size(x(i).df,2)
        for z = 1:size(y(i).df,2)
            [R P] = corrcoef(x(i).df(:,j),y(i).df(:,z));
            w(i).r(z,j) = R(1,2);
            w(i).p(z,j) = P(1,2);
        end
    end
end