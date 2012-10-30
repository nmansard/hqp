function [h,Y] = up(kup,cup,bound,h,Y,THR);

% up    actives one constraint and modifies consequently the HCOD.
%% Synopsis:
%    [h Y] = up(kup,cup,bound,h,Y)
%    [h Y] = up(kup,cup,bound,h,Y,THR)
%    [h Y] = up(kup,cup,"right",h,Y)
%    [h Y] = up(kup,cup,"left",h,Y)
%% Input:
%    kup   level where the constraint is going to be activated
%    cup   id of the constraint of level kup to be activated
%    bound (among 1="left" for upper bound and 2="right" for lower bound).
%               Bound to be activated.
%    h     "h" structure storing all the HQP data.
%    Y     right basis of the HCOD.
%    THR   epsilon threshold used to decide the rank-deficiency of a matrix. 
%% Output:
%    h,Y   the function modifies the input h and Y and returns them.
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- and Adrien Escande -- JRL/CNRS -- cf. COPYING.LESSER
%

% --- DEFAULT ARGUMENTS --------------------------------------------------------
nin = nargin
if nin==5
    % Default argument for the HCOD threshold.
    THR=1e-8;
    nin=nin+1;
end
% ---------------------------------------------------------------------

addpath('utils');
p = length(h);
nh = size(Y,1);

if isa(bound,'char')
    if strcmp(bound,'left') bound=1;
    elseif strcmp(bound,'right') bound=2;
    else print('Error with the bound arg.'); bound =1;
    end
end

% --- Current stage ------------------------------------------------------------
hk=h(kup);
iw=hk.iw; im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;

% Add a new row to H and a new dimension in W.
hk.iw                   = [iw hk.fw(1)];
iw                      = [iw hk.fw(1)];
hk.fw                   = hk.fw(2:end);
hk.im                   = [im hk.fm(1)];
im                      = [im hk.fm(1)];
hk.fm                   = hk.fm(2:end);
hk.H(im(end),:)         = hk.A(cup,:)*Y;
hk.W(iw(end),:)         = 0;
hk.W(:,im(end))         = 0;
hk.W(iw(end),im(end))   = 1; 
hk.m                    = m+1;
m                       = m+1;
hk.active(end+1,1)      = cup;
hk.activeb(end+1,1)     = cup + (bound==2)*hk.mmax;
hk.bound(cup)           = bound;

% 1.a Compute the rank of the new row, Eq (81).
%   Find the first element of the tail of the row such that the norm of the
%   tail is higher than the threshold: rup = max{ r, norm(ML(end,r:end))>THR }.
rup = nh+1-find(  cumsum(flipdim( hk.H(im(end),:).^2, 2 )) > THR^2,  1  );

% 1.b modify the decomposition of level k.
if rup<=ra
    % 1.b.TRUE (case V.B.2) The new row does not increase the rank of the level.
    if rup>rp
        % There is nonzero elements below L, Eq (82).
        for i=rup-rp:-1:1
            Wi          = givens( hk.H(im(:),rp+i),n+i,m )';
            hk.H(im,:)  = Wi*hk.H(im,:);
            hk.W(iw,im) = hk.W(iw,im)*Wi';
        end
        % else: There is only zeros below L, Eq (83): nothing to do.
    end
    % Switch H row.
    hk.im = [im(end) im(1:end-1)];
    hk.n  = n+1;
    h(kup)=hk; clear hk;
    return;
end

% 1.b.FALSE (case V.B.3) The new row does not increase the rank of the level.
%     Nullify the tail of the row: 
Yup = eye(nh);
for i=rup-1:-1:ra+1  % The loop corresponds to Eq (79).
    Yi               = givens(hk.H(im(end),:),i,i+1);
    hk.H(im(end),:)  = hk.H(im(end),:)*Yi;
    Yup              = Yup*Yi;
end
hk.r  = r+1;
hk.ra = ra+1;

h(kup)=hk; clear hk;
clear iw im r n ra rp m;


% --- Propagation --------------------------------------------------------------
% If Yup exists, apply it to the above levels (sec. V.B.4)
for k=kup+1:p
    massert( not(exist('hk')) ,'Error, hk should have been cleared.');
    hk=h(k);
    iw=hk.iw; im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;

    hk.H(im,:) = hk.H(im,:)*Yup;
    
    if rup<rp+1
        % (case V.B.7): Rank already lost, nothing to do.
    elseif rup>ra
        % (Case V.B.5): Rank lost later, L simply shifts on the right.
        hk.rp = rp+1;
        hk.ra = ra+1;
    else
        % (case V.B.6) Aup is colinear to one row of this level: rank lost here.
        rdef = rup-rp; % Row to be removed from L.
        for i=rdef:-1:2 % The loop corresponds to Eq (79).
            Wi             = givens(hk.H(im(:),rp+i),n+i-1,n+rdef)';
            hk.H(im(:),:)  = Wi*hk.H(im(:),:);
            hk.W(iw,im)    = hk.W(iw,im)*Wi';
        end
        
        hk.rp = rp+1;
        hk.r  = r-1;
        hk.n  = n+1;
        hk.im = im( [ n+rdef  1:n n+1:n+rdef-1 n+rdef+1:m] ); 
    end
    
    h(k)=hk; clear hk;
    clear iw im r n ra rp m;
end
Y=Y*Yup;
