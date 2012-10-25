function L = act_L(h,k);

if nargin==2
    if strcmp(k,'all')
        L={};
        for k=1:length(h)
           L{k} = act_L(h,k); 
        end
      
    else
        ir = h(k).im(h(k).n+1:end);
        L = h(k).H(ir,h(k).rp+1:h(k).ra);
    end
else
    L=[];
    for k=1:length(h)
        Lk = act_L(h,k); rk=size(Lk,1);
        
        L(end+1:end+rk,end+1:end+rk)= Lk;
    end
end
