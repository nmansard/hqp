function res = dispaset(h)

% Return a string displaying the active set in h.

p=length(h);
res='';
for k=1:p
    hk=h(k);
    str=' => ';
    for c=sort(hk.active')
        if hk.btype(c)==1 sgn='=';
        elseif hk.bound(c)==1 sgn='-';
        elseif hk.bound(c)==2 sgn='+';
        end
        str=sprintf('%s %s%d',str,sgn,c);
    end
    res=sprintf('%s\n%s',res,str);
end
