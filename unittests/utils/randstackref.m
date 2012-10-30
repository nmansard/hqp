function [ nh p m r seed ] = randstackref(seed)

% randstackref randomly generates the size of a HQP problem. These size
% parameters can then be passed to randstack to get a random HQP.
% 
%% Synopsis:
%  [ nh p m r seed ] = randstack()
%  [ nh p m r ]      = randstack(seed)
%
% Input:
%  seed   if given, the seed is used to initialize the random generator. 
%            Otherwise, a seed is randomly chosen and returned.
%
%% Output:
%   nh       is the size of the parameter space (number of columns of A).
%   p        is the number of levels.
%   m        is the size (number of rows) of each levels.
%   r        is the rank of each levels (rank of Ak Z_{k-1} when all the
%               constraints are active).
%   seed     is the randomly chosen seed if it was not given in input.
%
%

% --------------------------------------------------------------------
if nargin==0
    seed = floor(rand*10000);
end
% --------------------------------------------------------------------

rand_init(seed);

nh=floor(exprnd(1)*20+5);
maverage = max(nh * (randn*0.2+0.25),4);
p=max(1,floor(nh/maverage));
m=[];
r=[];
m = ceil(max(0.2,randn(p,1)*0.5+1)*maverage);
r = floor(max(0.2,randn(p,1)*0.5+0.8)*maverage);
r = min(m,r); r=max(r,1);
if sum(r)>nh
    r=floor(r/sum(r)*nh);
end
