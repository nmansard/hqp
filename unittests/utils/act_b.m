function b = act_b(h,k);

if nargin==2
    if strcmp(k,'all')
        b={};
        for k=1:length(h)
            b{k} = act_b(h,k);
        end
      
    else
        ia = h(k).active;
        b = h(k).b(h(k).activeb);
    end
else
    b=[];
    for k=1:length(h)
        ia = h(k).active;
        b = [ b; act_b(h,k)];
    end
end
