function [Fr] = pos2indx(pos,nms,ms)





upos = unique(pos);

for w = 1:size(upos,1)
    switch upos(w)
        case 0
            Fr(1,1) = mean(nms(pos == upos(w))); % NMS
            Fr(1,2) = mean(ms(pos == upos(w))); % MS
        case 1
            Fr(2,1) = mean(nms(pos == upos(w))); % NMS
            Fr(2,2) = mean(ms(pos == upos(w))); % MS
        case 2
            Fr(3,1) = mean(nms(pos == upos(w))); % NMS
            Fr(3,2) = mean(ms(pos == upos(w))); % MS
        case 3
            Fr(4,1) = mean(nms(pos == upos(w))); % NMS
            Fr(4,2) = mean(ms(pos == upos(w))); % MS
        case 4
            Fr(5,1) = mean(nms(pos == upos(w))); % NMS
            Fr(5,2) = mean(ms(pos == upos(w))); % MS
        case 5
            Fr(6,1) = mean(nms(pos == upos(w))); % NMS
            Fr(6,2) = mean(ms(pos == upos(w))); % MS
    end
    
end