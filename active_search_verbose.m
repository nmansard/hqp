function [ primal dual h Y xtrack ] = active_search_verbose(A,b,btype,aset_init, aset_bound,THR);

% active_search  solves a iHQP given in input using a hiearchical active search
%                   method and prints all the step of the computation. This 
%                   function corresponds to Alg. 5. It is strictly similar to 
%                   active_search.m, except for the traces.
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
%    aset_init     if specify, is the initial active search (cell of vector
%                    of ranges of b rows)
%    aset_bounds   are the corresponding bounds (=,<=,>=,<=<=). 
%% Output:
%    primal        primal optimum answering to the HQP.
%    dual          corresponding dual optimum.
%    h             "h" structure storing optimal active set and the
%                      corresponding HCOD.
%    Y             right basis of the HCOD "h".
%    THR           is the threshold used to test the positivity.
%
% The function corresponds to Alg. 5. To avoid an extra computation of the
% primal at each iteration of the outer loop of the paper, it is performed
% inside the "for k=kcheck+1:p+1" loop, and the current level of the outer
% loop is stored inside "kcheck".
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- cf. COPYING.LESSER
%

addpath('utils');
p=length(A);
nh=size(A{1},2);

% --- DEFAULT ARGUMENTS --------------------------------------------------------
nin = nargin
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

VERBOSE = true

% --- Initial HCOD
% The active set and the HCOD are stored in the cell "h". See the hcod
% documentation for details.
[h Y]   = hcod(A,b,btype,aset_init,aset_bound);
y0 = zeros(nh,1);
x0 = zeros(nh,1);
kcheck  = 0;               % level of the "outer" loops whose multiplier
                           % have already been computed and tested.
iter    = 1;               % Number of iterations of the active search.
dual    = {};              % Stored Lagrangian multipliers when kcheck grows.
xtrack  = {};

while kcheck<=p
    [ y1 x1 ]              = ehqp_primal(h,Y);              % Alg 5#8
    xtrack{1,iter} = x0;
    xtrack{2,iter} = x1;
    [viol tau cst]         = step_length(x0,x1,h,Y);        % Alg 5#10
    if viol
        if VERBOSE dispcst('Violation',iter,cst); end
        x0                 = (1-tau)*x0+tau*x1;             % Alg 5#11
        [h Y]              = up(cst(1),cst(2),cst(3),h,Y);  % Alg 5#14
    else
        x0=x1;
        for k=kcheck+1:p+1
            lambda         = ehqp_dual(k,y1,h,Y);           % Alg 5#18
            [need cst ]    = check_mult(lambda,h,Y);        % Alg 5#20
            if need                                         % Alg 5#21
                if VERBOSE 
                    level=cst(1);
                    c_id=h(level).active(cst(2));
                    c_bound=h(level).bound(c_id);
                    dispcst('Suboptimal',iter,[level c_id c_bound]);
                end
                [h Y]      = down(cst(1),cst(2),h,Y);       % Alg 5#22
                break;
            end

            h              = freeze(lambda,h);              % Alg 5#25
            dual{k}        = lambda;
            kcheck         = k;
            if VERBOSE disp(sprintf('Validation of level %d',kcheck)); end
        end
    end
    iter=iter+1;
end

primal = x1;
xtrack{1,end+1} = x1;
xtrack{2,end} = x1;
