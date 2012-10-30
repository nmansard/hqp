function [ A b btype nh p m r Au bu  ] = randstack( nh,p,m,r,svbound );

% randstack randomly generates a hierarchical problem from the size of each 
% level and the bounds of the singular values.
% 
%% Synopsis:
%  [ A b btype nh p m r Au bu ] = randstack()
%  [ A b btype nh p m r Au bu ] = randstack(nh)
%  [ A b btype nh p m r Au bu ] = randstack(nh,p,m,r)
%  [ A b btype nh p m r Au bu ] = randstack(nh,p,m,r,svbound)
%
%% Input:
%   nh       is the size of the parameter space (number of columns of A).
%   p        is the number of levels.
%   m        is the size (number of rows) of each level.
%   r        is the rank of each level (rank of Ak Z_{k-1} when all the
%               constraints are active).
%   svbound  is the bound of the nonzero singular values of Ak Z_{k-1}
%               when all the constraints are active.
%
%% Output:
%   A,b        are the square-root Hessian and gradient of each level of
%                 constraints.
%   btype      is the cell of the type of constraints (=,<=,>=,<=<=).
%   nh,p,m,r   are the problem size given in input.
%   Au,bu      are a copy of A and b under a matrix form (instead of a
%                 cell).
%
%
    
% --------------------------------------------------------------------
nin = nargin;
if nin == 0
    nin=1;
    nh=12;
end

if nin == 1
    if nh==12
        % Size of the problem
        nh=12;
        % Number of levels
        p=4;
        % Size of each level
        m=[3 4 3 4];
        % Rank of each level
        r=[ 2 3 2 3 ];
    elseif m==8
        % Number of levels
        p=2;
        % Size of each level
        m=[3 4];
        % Rank of each level
        r=[ 2 3 ];
    end

    nin = 4;
end

if nin==4
    svbound = [0.5 1.5];
end
% --------------------------------------------------------------------

Au=[]; bu=[];
for k=1:p
    rak = sum(r(1:k-1));
    mk=m(k); rk=r(k); nk=mk-rk; zk = nh-rak-rk;
    
    % Generate the bounds of the problem, as b = [ x1 x2 ] with x1<x2.
    b{k}=sort(((-.5+rand(mk,2))*2),2);
    % Generate the type of the constraints: nEq equalities first, and then
    % bounds among <=,>=,<=<= .
    nEq=round(rand*m(k));
    btype{k} = [ ones(nEq,1) ; ceil(rand(m(k)-nEq,1)*3)+1 ];
   
    % Generate a matrix of the form W * [ N 0 0 ; M S 0 ] with S a nonzero
    % diagonal a W full rank.
    Sk = diag(rand(rk,1)*diff(svbound)+svbound(1));
    Ak = [ rand(mk,rak) [ Sk;zeros(nk,rk) ] zeros(mk,zk) ];
    Wk = mrand(mk,mk,[1 1]);
    Au = [Au; Wk*Ak];
end

% Au is finally of the form W*H*Y, with W a block diagonal of full rank
% matrices, Y a full rank matrix, and H as diaongal blocks instead of the Lk
% of the HCOD.
Au = Au*mrand(nh,nh,[1 1]);

% Divide the Au matrix into a cell of p matrices A_1 ... A_p.
for k=1:p
    maprec = sum(m(1:k-1)) ;
    A{k} = Au( maprec+1:maprec+m(k) , : );
    bu=[bu;b{k}];
end

end


% --------------------------------------------------------------------
% --------------------------------------------------------------------
% --------------------------------------------------------------------
function [M] = mrand(n,m, sbound)

  % Generate a matrix whose svd is inside sbound (default is [.5 1.5]).

  if nargin==2
      sbound = [0.5 1.5];
  end

  smin=sbound(1);
  slength = diff(sbound);

  r=min(n,m);

  [U S V ] = svd( rand(n,m));
  S(1:r,1:r) = diag(rand( min(n,m),1 )*slength+smin );

  M=U*S*V';

end