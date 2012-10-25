function [R] = givens(x,i,j)

% Givens    compute the Givens rotation R that nullify x(j) using x(i):
%              (R*x)(j) == (x'*R')(j) == 0
%% Synopsis:
%      R = givens(x,i,j)
%% Input:
%      x    vector used to compute the angle of the rotation
%      i    index of the coefficient of x used as a pivot
%      j    index of the coefficient of x to be nullified by R.
%% Output:
%      R    the Givens rotation G[i,j,theta] with theta = atan(xj,xi)
%
%

m=length(x);
R=eye(m);
xi=x(i);xj=x(j);

if abs(xj)<1e-9
    c=1; s=0;
elseif abs(xj)>abs(xi)
    theta=-xi/xj; s=1/sqrt(1+theta*theta); c=s*theta;
else
    theta=-xj/xi; c=1/sqrt(1+theta*theta); s=c*theta;
end

R(i,i) =  c;
R(i,j) = -s;
R(j,j) =  c;
R(j,i) =  s;
