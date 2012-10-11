function [aset b] = dummybound(aset);
% return a cell of array full of ones that have the same size than the arg

p=length(aset);
if strcmp(class(aset),'struct')
    h=aset; aset={};
    m=[h.mmax];
end
if strcmp(class(aset),'double')
    m=aset; aset={};
end

for k=1:p
    if length(aset)<k aset{k} = 1:m(k); end
    if size(aset{k},1)==1 aset{k}=aset{k}'; end
    b{k} = ones(length(aset{k}),1);
end

