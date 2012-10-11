function A = _A(h,k);

if nargin==2
    if strcmp(k,'all')
        A={};
        for k=1:length(h)
            ia = h(k).active;
            A{k} = h(k).A(ia,:);
        end
      
    else
        ia = h(k).active;
        A = h(k).A(ia,:);
    end
else
    A=[];
    for k=1:length(h)
        ia = h(k).active;
        A = [A;h(k).A(ia,:)];
    end
end
