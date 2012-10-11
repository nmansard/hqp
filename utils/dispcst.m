function str = dispcst(txt,iter,cst);
% return the description of the current constraint

str=sprintf('%s(%d):',txt,iter);

for k=1:size(cst,1)
    if cst(k,3)==1            sgn = '-';
    elseif cst(k,3)==2        sgn = '+';
    else                      sgn = '' ;
    end

    str = sprintf('%s %d:%s%d',str,cst(k,1),sgn,cst(k,2));
end

if nargout==0
    disp(str);
end
