function H = _H(h,k);

if nargin==2
    if strcmp(k,'all')
        H={};
        for k=1:length(h)
           H{k} = _H(h,k); 
        end
      
    else
        im = h(k).im;
        H = h(k).H(im,1:h(k).ra);
    end
else
    H=[];
    for k=1:length(h)
        Hk = _H(h,k);
        if k>1 & (size(H,2)<size(Hk,2)) 
            H(end,size(Hk,2)) = 0;
        end
        H= [H;Hk];
    end
end
