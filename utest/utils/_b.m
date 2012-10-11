function b = _b(h,k);

if nargin==2
    if strcmp(k,'all')
        b={};
        for k=1:length(h)
            b{k} = _b(h,k);
        end
      
    else
        ia = h(k).active;
        b = h(k).b(h(k).activeb);
    end
else
    b=[];
    for k=1:length(h)
        ia = h(k).active;
        b = [ b; _b(h,k)];
    end
end
