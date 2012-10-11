function slack = hdist(x,h);

% compute the slack of all active constraints wrt solution x. Check btw that
% the unactive constraints are verified.
% For each level, return an array of 4 rows and mmax columns. Each component
% of the level is informed by one column:
%    * Row 1: 0 if inactive, otherwise, 1 if twin, 3 if inf, 4 if double.
%    * Row 2: distance to the upper bound
%    * Row 3: distance to the lower bound
%    * Row 4: type of the constraint, if active.

p=length(h);
Etwin = 1;Edouble = 2;Einf = 3;Esup = 4;Enone=0;
EPS=1e-3;
slack={};
for k=1:p
    hk=h(k);
    Ax=hk.A*x;
    slack{k}=[];
    for c=1:hk.mmax
        slack{k}(4,c) = hk.btype(c);

        if ismember(c,hk.active)

            slack{k}(1,c) = hk.bound(c);
            
            if     hk.btype(c)==Etwin
                slack{k}(2:3,c)=hk.b(c,1)-Ax(c);
            elseif hk.btype(c)==Einf
                slack{k}(2,c)=Ax(c)-hk.b(c,1);
                slack{k}(3,c)=0;
            elseif hk.btype(c)==Esup
                slack{k}(2,c)=0;
                slack{k}(3,c)=hk.b(c,2)-Ax(c);
            elseif hk.btype(c)==Edouble
                slack{k}(2,c)=Ax(c)-hk.b(c,1);
                slack{k}(3,c)=hk.b(c,2)-Ax(c);
            end                
            
        else
            slack{k}(1,c) = 0;
            
            if     hk.btype(c)==Etwin
                slack{k}(2:3,c)=nan;
            elseif hk.btype(c)==Einf
                slack{k}(2,c)=Ax(c)-hk.b(c,1);
                slack{k}(3,c)=0;
            elseif hk.btype(c)==Esup
                slack{k}(2,c)=0;
                slack{k}(3,c)=hk.b(c,1)-Ax(c);
            elseif hk.btype(c)==Edouble
                slack{k}(2,c)=Ax(c)-hk.b(c,1);
                slack{k}(3,c)=hk.b(c,2)-Ax(c);
            end                
        end
    end

end


