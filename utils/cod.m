function [ W L Q rankA ] = cod(A,THR,smallL);

% cod   compute the lower-triangular complete orthogonal decomposition. The
% epsilons smaller than the input threshold are neglected.
% a
%% Synopsis:
%     [ W L Q ]      = cod(A)
%     [ W L Q ]      = cod(A,THR)
%     [ W L Q rank ] = cod(A)
%     [ W L Q rank ] = cod(A,THR)
%     [ W L Q ]      = cod(A,THR,1)
%
%% Input:
%    A        is the matrix to decompose
%    THR      is the threshold below which rows are neglected.
%    smallL   if specified to 1, the non zero part of L is returned (in which 
%               case the product W*L*Q is not consistent).
%% Output:
%    W,L,Q    returns the COD decomposition W,L,Q such that W*L*Q==A, W and
%               Q being orthonormal matrices.
%    rankA    rank of A wrt "THR", ie size of L.
%
%

% --- DEFAULT ARGUMENTS --------------------------------------------------------
if nargin<2
    THR=1e-3;
end
% ---------------------------------------------------------------------

% 1. Perform a column-pivoting QR decomposition of A' (ie a row-pivoting LQ
% decomposition of A). Eq (77).
[Q R E]=qr(A');
L=R'; %            A == E L Q'

% 2. Compute the rank of A.
rankA=0;
for i=1:min(size(A))
  if(abs(L(i,i))>THR) rankA=rankA+1;
  end;
end;
[m n]=size(A);

% 3. Compute the left basis of the COD.
W=eye(m);

% This loop corresponds to eq. (II-15).
if rankA<m
  for j=rankA:-1:1
    for i=m:-1:rankA+1
        U=givens(L(:,j),j,i);
        L=U*L;  W=W*U';
    end
  end
end

% Permute so as to have L=[0;L0] instead of L=[L0;0].
L=[L(rankA+1:end,:); L(1:rankA,:)];
W=E*[W(:,rankA+1:end) W(:,1:rankA)];

% Shorten L to L0 if the second argument cod(A,0) is here.
if nargin==3
    nullA=m-rankA;
    L=L(nullA+1:end,1:rankA);
end

