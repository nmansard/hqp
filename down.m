function [h,Y] = down(kdown,rdown,h,Y,THR);

% down    deactivates one constraint and modifies consequently the HCOD.
%% Synopsis:
%    [h Y] = down(kdown,rdown,h,Y)
%    [h Y] = down(kdown,rdown,h,Y,THR)
%% Input:
%    kdown   level where the constraint is going to be deactivated
%    rdown   row number of the constraint of level kdown to be deactivated
%    h       "h" structure storing all the HQP data.
%    Y       right basis of the HCOD.
%    THR      epsilon threshold used to decide the rank-deficiency of a matrix.
%% Output:
%    h,Y   the function modifies the input h and Y and returns them.
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- cf. COPYING.LESSER
%

% --- DEFAULT ARGUMENTS --------------------------------------------------------
nin = nargin
if nin==4
    % Default argument for the HCOD threshold.
    THR=1e-8;
    nin=nin+1;
end
% ---------------------------------------------------------------------

addpath('utils');
p = length(h);
nh = size(Y,2);

% ------------------------------------------------------------------------------
% 1. Remove the line from current stage (Sec V.C.1)
hk=h(kdown);
iw=hk.iw; im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;

% 1.a Set a "1" in the <ldown> row of W: Eq (89).
% W is transformed to:
%   W = [ W_a  0  W_b ]
%       |  0   1   0  |
%       [ W_c  0  W_d ]

% Search for the first non-zero of W(row,:).
pivot = find(abs(hk.W(iw(rdown),im))>THR,1,'first');

% Uses Givens rotations to nullify the tail of the row.
for i=pivot+1:m
    if abs(abs(hk.W(iw(rdown),im(pivot)))-1)<THR^2/2
        break;
    end
    Wdown         = givens(hk.W(iw(rdown),im),pivot,i);
    hk.W(iw,im)   = hk.W(iw,im)*Wdown;
    hk.H(im,1:ra) = Wdown'*hk.H(im,1:ra);
end

% 1.b remove the corresponding row.
hk.fw      = [iw(rdown) hk.fw];
hk.fm      = [im(pivot) hk.fm];
hk.active  = hk.active([1:rdown-1 rdown+1:m]);
hk.activeb = hk.activeb([1:rdown-1 rdown+1:m]);
hk.iw      = iw([1:rdown-1 rdown+1:m]);
iw         = iw([1:rdown-1 rdown+1:m]);
hk.im      = im([1:pivot-1 pivot+1:m]);
im         = im([1:pivot-1 pivot+1:m]);
hk.m       = m-1;
m          = m-1;

% Bug in octave when a vector becomes []: patch below.
if length(hk.active)==0
    hk.active=zeros(0,1);
    hk.activeb=zeros(0,1);
end

% 1.c. CASE 1: no modification rank to be reported.
if pivot<=hk.n
    hk.n = n-1;
    n    = n-1;
    h(kdown)=hk; clear hk;
    return;
end

% 1.c CASE 2: L is now upper Hessenberg and should be corrected.
Ydown=eye(nh);
for i=1:hk.r-1
    R=givens(hk.H(im(n+i),:),rp+i,rp+i+1);
    Ydown=Ydown*R;
    hk.H(im,:)=hk.H(im,:)*R;
end

r     = hk.r-1;
hk.r  = hk.r-1;
ra    = hk.ra-1;
hk.ra = hk.ra-1;
h(kdown)=hk; clear hk;

% ------------------------------------------------------------------------------
% 2. Propagation of Ydown to the levels k+1 .. p (Sec V.C.1)

% noMoreRankChange will be set to 1 when one level increases its rank (last
% paragraph of Sec. V.C.2).
noMoreRankChange=0;

for k=kdown+1:p
    massert(not(exist('hk')),'Error, hk should have been cleared.');
    hk = h(k);
    iw=hk.iw; im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;

    % 2.a Apply Ydown
    hk.H(im,:) = hk.H(im,:)*Ydown;

    if( noMoreRankChange )
        h(k)=hk; clear hk;
        continue;
    end
    
    % Check among the m1..mn of Eq (90) if one is nonzero.
    if norm(hk.H(im(1:n),hk.rp))>THR 

        % 2.b.CASE 1: one of the m1..mn is nonzero
        
        % Find the first nonzero m_i.
        prom = find( abs(hk.H(im(1:n),rp)),1 );

        % Nullify all the m_i from m_prom.
        for i=prom+1:hk.n
            Wdown         = givens(hk.H(im,rp),prom,i)';
            hk.H(im,1:ra) = Wdown*hk.H(im,1:ra);
            hk.W(iw,im)   = hk.W(iw,im)*Wdown';
        end

        % Reorganize H_k.
        hk.im             = im([1:prom-1 prom+1:n prom n+1:m]);
        im                = im([1:prom-1 prom+1:n prom n+1:m]);
        hk.rp             = rp-1;
        rp                = rp-1;
        hk.r              = r+1;
        r                 = r+1;
        hk.n              = n-1;
        n                 = n-1; 
        noMoreRankChange  = 1;
    else
        % 2.b.CASE 2: all the m_i are zero, L_k is upper hessenberg.
        for i=1:r
            R             = givens(hk.H(im(n+i),:),rp-1+i,rp-1+i+1);
            Ydown         = Ydown*R;
            hk.H(im,:)    = hk.H(im,:)*R;
        end
        hk.rp = rp-1;
        rp    = rp-1;
        hk.ra = ra-1;
        ra    = ra-1;
    end
    
    h(k)=hk; clear hk;
end

Y=Y*Ydown;
