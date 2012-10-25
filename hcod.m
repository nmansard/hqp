function [h Y] = hcod(A,b,btype,active,bounds,EPS);

% hcod      computes the HCOD decomposition of A, and initialize the structure
%              h storing the problem data.
%% Synopsis:
%    [h Y] = hcod(A,b,btype)
%    [h Y] = hcod(A,b,btype, active,bounds )
%    [h Y] = hcod(A,b,btype, active,bounds,EPS)
%
%% Input:
%    A,b      square-root Hessian and gradient of the HQP problem.
%    btype    bound types (=,<=,>=,<=<=) of the constraints
%    active   active set for which the HCOD should be computed.
%    bounds   active bounds of the active set
%    EPS      epsilon threshold used to decide the rank-deficiency of a matrix. 
%    
%% Output:
%    h        h structure storing the whole problem, the active set and the
%                decomposition. h is detailed below.
%    Y        the right basis (common to all the levels and then not stored in 
%                the cell h).
%
%% Initialization of the h structure.
% === PROBLEM DEFINITION ===
% h{k}.A         is a copy of the A_k matrix.
% h{k}.b         is a copy of the b_k vectors.
% h{k}.btype     is the types (=,<=,>=,<=<=) of the constraints A,b.
%
% === INDICES ===
% h{k}.mmax      is the number of constraints = rows(A_k).
% h{k}.m         is the number of active constraints.
% h{k}.r         is the rank of the level = rank(A_k Z_{k-1}).
% h{k}.n         is the number of rank-deficient constraints = m-r.
% h{k}.ra        is the number of columns of H_k = sum(i=1:k) r_i
% h{k}.rp        is the number of columns of M_k = sum(i=1:k-1) r_i
%
% === HCOD MATRICES ===
% Arrays are allocated for W and H with maximum possible size. They are then
% accessed with the corresponding index vector.
% h{k}.iw        is the index vector of the rows of W.
% h{k}.im        is the index vector of H rows and W columns.
% h{k}.W         is the memory where W_k is stored.
%                  * W_k = h{k}.W(iw,im), with im=h{k}.im and iw=h{k}.iw
%                  * V_k = h{k}.W(iw,im0), with im0=im(1:n) and n=h{k}.n
%                  * U_k = h{k}.W(iw,im1), with im1=im(n+1:end)
% h{k}.H         is the memory where H_k is stored.
%                  * H_k = h{k}.H(im,1:ra), with ra=h{k}.ra
%                  * N_k = h{k}.H(im0,1:rp), with r=h{k}.rp
%                  * M_k = h{k}.H(im1,1:rp)
%                  * L_k = h{k}.H(im1,rp+1:ra)
% h{k}.fw        is the vector of free rows of W.
% h{k}.fh        is the vector of free rows of H.
%
% === ACTIVE SET ===
% h{k}.active    is the map of the active constraints in W H Y':
%                  * (W*H*Y')(i,:) = A(active(i),:)
% h{k}.bound     is the type of the constraints (upper or lower bounds, when relevant).
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- cf. COPYING.LESSER
%

addpath('utils');
p  = length(A);    % Number of level
nh = size(A{1},2); % Parameter size

% --- DEFAULT ARGUMENTS --------------------------------------------------------
nin = nargin
if nin==2
    % Default argument for "btype": all the constraints are equalities.
    for k=1:p
        btype{k}=1:size(A{k},1);
    end
    nin=nin+1;
end

if nin==3
    % Default argument for "active": all the constraint are active.
    for k=1:p
        active{k}=1:size(A{k},1);
    end
    nin=nin+1;
end

if nin==4
    % Default argument for "bounds": all the active constraints are the lower
    % bounds.
    for k=1:p
        bounds{k} = ones(length(active{k}),1);
    end
    nin=nin+1;
end

if nin==5
    % Default argument for the HCOD threshold.
    EPS=1e-8;
    nin=nin+1;
end
% ---------------------------------------------------------------------

for k=1:p

    if size(active{k},1)==1 active{k}=active{k}'; end
    if size(bounds{k},1)==1 bounds{k}=bounds{k}'; end

    hk.A=A{k};
    hk.b=b{k};
    hk.btype=btype{k};
    hk.mmax = size(hk.A,1);

    hk.active = active{k};
    hk.bound = zeros(hk.mmax,1);
    hk.bound(hk.active) = bounds{k};
    hk.activeb = active{k} + (bounds{k}==2)*hk.mmax;
    hk.freeze = zeros(hk.mmax,1);
    
    hk.W = eye(hk.mmax);
    hk.H = zeros(hk.mmax,nh);
    hk.m = length(hk.active);

    hk.iw = 1:hk.m;
    hk.im = 1:hk.m;
    % W = hk.W(iw,im) -- NM = hk.H(im,1:rp) -- L = hk.H(im(n+1:s),rp+1:ra)
    hk.fw = hk.m+1:hk.mmax;
    hk.fm = hk.m+1:hk.mmax;

    hk.n = 0;
    hk.r = 0;
    hk.rp= 0;
    hk.ra= 0;
    
    h(k) = hk;
end

% --- HCOD LOOP -------------------------------------------------------

% Decompose the first level using a simple COD.
k=1;
hk=h(k);
if hk.m>0
    [ W L Y hk.r]     = cod(A{k}(active{k},:),EPS);
    hk.W(hk.iw,hk.im) = W;
    hk.H(hk.im,:)     = L;
    hk.ra             = hk.r;
    hk.n              = hk.m-hk.r;
else
    hk.ra             = 0;
    Y                 = eye(nh);
end
h(k)=hk;

% Compute the HCOD for the other levels.
for k=2:p
    hk=h(k);
    hk.rp=h(k-1).ra;
    m=hk.m; im=hk.im; iw=hk.iw; rp=hk.rp;

    if m>0
        AY = hk.A(hk.active,:)*Y(:,rp+1:end);

        % Compute W, L and Yk.
        [W H Yk hk.r]     = cod(AY,EPS);
        hk.W(iw,im)       = W;
        hk.H(im,rp+1:end) = H;
        hk.ra             = rp+hk.r;
        hk.n              = hk.m-hk.r;

        % Compute M and N.
        hk.H(hk.im,1:hk.rp) = hk.W(hk.iw,hk.im)'*hk.A(hk.active,:)*Y(:,1:hk.rp);

        % Apply Yk on Y.
        Y = [ Y(:,1:rp) Y(:,rp+1:end)*Yk ];
    else
        hk.ra = rp;
    end
    h(k)=hk;
end

