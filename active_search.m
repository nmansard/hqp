function [ primal dual h Y ] = active_search(A,b,btype,aset_init, aset_bound,THR);

% active_search  solves a iHQP given in input using a hiearchical active search
%                   method. This function corresponds to Alg. 5.
%% Synopsis:
%     primal              = active_search(A,b,btype);
%     primal              = active_search(A,b,btype,aset_init,aset_bound);
%     [ primal dual ]     = active_search(A,b,btype,aset_init,aset_bound);
%     [ primal dual h Y ] = active_search(A,b,btype,aset_init,aset_bound);
% 
%% Input:
%    A             is <p> cells of <m_k>x<nh>  matrix .
%    b             is the <p> cells of <m_k>x2 matrices (first column is dummy 
%                    for sup bounds, second column is dummy for inf and twin 
%                    bounds). 
%    btype         are the type (=,<=,>=,<=<=) for each b rows.
%    aset_init     if specify, is the initial active search (cell of vector of 
%                    ranges of b rows)
%    aset_bounds   are the corresponding bounds (=,<=,>=,<=<=). 
%% Output:
%    primal        primal optimum answering to the HQP.
%    dual          corresponding dual optimum.
%    h             "h" structure storing optimal active set and the
%                      corresponding HCOD.
%    Y             right basis of the HCOD "h".
%    THR           is the threshold used to test the positivity.
%
% The function corresponds to Alg. 5 of the paper Hierarchical Quadratic
% Programming (Part I). To avoid an extra computation of the primal at each
% iteration of the outer loop of the paper, it is performed inside the "for
% k=kcheck+1:p+1" loop, and the current level of the outer loop is stored
% inside "kcheck".
%
% Copyright Nicolas Mansard -- LAAS/CNRS
%    -- and Adrien Escande -- JRL/CNRS
%    -- cf. LICENSE.txt
%

addpath('utils');
p=length(A);
nh=size(A{1},2);

% --- DEFAULT ARGUMENTS --------------------------------------------------------
nin = nargin;
if nin==3
    % The initial active set only contains equality constraints.
    [ aset_init aset_bound ] = initset(btype);
    nin=5;
else
    % Check the active set and activate equality constraints.
    [ aset_init aset_bound ] = initset(btype,aset_init,aset_bound);
end

if nin==5
    THR = 1e-8;
end
% ---------------------------------------------------------------------

% --- Initial HCOD
% The active set and the HCOD are stored in the cell "h". See the hcod
% documentation for details.
[h Y]   = hcod(A,b,btype,aset_init,aset_bound);
y0 = zeros(nh,1);
x0 = zeros(nh,1);
kcheck  = 0;               % level of the "outer" loops whose multiplier
                           % have already been computed and tested.
iter    = 0;               % Number of iterations of the active search.
dual    = {};              % Stored Lagrangian multipliers when kcheck grows.

while kcheck<=p
    [ y1 x1 ]              = ehqp_primal(h,Y);              % Alg 5#8
    [viol tau cst]         = step_length(x0,x1,h,Y);        % Alg 5#10
    if viol
        x0                 = (1-tau)*x0+tau*x1;             % Alg 5#11
        [h Y]              = up(cst(1),cst(2),cst(3),h,Y);  % Alg 5#14
    else
        x0=x1;
        for k=kcheck+1:p+1
            lambda         = ehqp_dual(k,y1,h,Y);           % Alg 5#18
            [need cst ]    = check_mult(lambda,h,Y);        % Alg 5#20
            if need                                         % Alg 5#21
                [h Y]      = down(cst(1),cst(2),h,Y);       % Alg 5#22
                break;
            end

            h              = freeze(lambda,h);              % Alg 5#25
            dual{k}        = lambda;
            kcheck         = k;
        end
    end
    iter=iter+1;
end

primal = x1;
