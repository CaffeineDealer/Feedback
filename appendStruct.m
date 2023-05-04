function [c] = appendStruct(s1,s2)


% Convert Struct to Cell
C1 = struct2cell(s1);
C2 = struct2cell(s2);
% Merg desired field
C = [C1{1,1}(:).xr C2{1,1}(:).xr];



end