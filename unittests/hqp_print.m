function str = hqp_print(h,filename);

% hqp_print prints a given HQP problem (with the "h" structure shape) into a
% file.
%% Synopsis:
%     str   = hqp_print(h)
%             hqp_print(h,filename)
% 
%% Input:
%   h         the h structure built by the hcod function.
%   filename  file where the problem should be written.
%% Output:
%   str       returns the string that may be written in the file.

p=length(h);
nh=size(h(1).A,2);
addpath('utils');
constants;

str=sprintf('variable size %d\n',nh);

for k=1:p
    hk=h(k);
    ie=find(hk.btype==Etwin)';
    ii=find(hk.btype!=Etwin)';
    ne=length(ie);
    ni=length(ii);

    str=sprintf('%slevel\nequalities %d\n',str,ne);
    for c=ie
        for i=1:nh
            str=sprintf('%s %.20f',str,hk.A(c,i));
        end
            str=sprintf('%s %.20f',str,hk.b(c,1));
   end
    
    
    str=sprintf('%s\ninequalities %d\n',str,ni);
    for c=ii
        bi=sprintf('%.20f',hk.b(c,1));        be=sprintf('%.20f',hk.b(c,2));
        
        if hk.btype(c)==Einf
            be='1e25';
        elseif hk.btype(c)==Esup
            bi='-1e25';
        end
        
        for i=1:nh
            str=sprintf('%s %.20f',str,hk.A(c,i));
        end
        str=sprintf('%s       %s %s\n',str,bi,be);
    end
    
end
str=sprintf('%send\n',str);

if nargin==2
    fileID = fopen(filename,'w');
    fprintf(fileID,'%s',str);
    fclose(fileID);
end
