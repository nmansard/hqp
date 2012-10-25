function [ need cst h  maxl ] = check_mult(lambda,h,Y,THR);

% check_mult searches a lexicographic positive multipliers that is not an
%               equality.
%% Synopsis:
%   [ need cst      ]   = check_mult(lambda,h,Y)
%   [ need cst maxl ]   = check_mult(lambda,h,Y)
%   [ need cst maxl h ] = check_mult(lambda,h,Y)
% 
%% Input:
%    lambda   multiplier to be tested.
%    h       "h" structure storing all the HQP data.
%    Y       right basis of the HCOD.
%    THR     is the threshold used to test the positivity.
%% Output:
%    need    returns "need = false" if no constraint satisfies the
%               lexicographic and bound-type properties. Otherwise, returns
%               "need = true" and the reference on the maximum.
%    cst     if need is true, return the reference to the constraint
%               corresponding of the maximum of the multiplier.
%    h       If need is false, then the function has a side effect: it
%               modifies the "freeze" field of the "h" structure. To account
%               for the side effect, h is returned.
%    maxl    the reached maximum corresponding to the cst constraint.
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- cf. COPYING.LESSER
%

% --- DEFAULT ARGUMENTS --------------------------------------------------------
if nargin<4
    THR=1e-8;
end
% ---------------------------------------------------------------------

kl=length(lambda);
nh=size(Y,2);
p=length(h);
constants;

maxl = 0; cst =[];

% --- POSITIVITY LOOP ------------------------------------------------
for k=1:kl
    massert( not(exist('hk')) ,'Error, hk should have been cleared.');
    hk=h(k);
    iw=hk.iw; im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;
    ia=hk.active;
    
    % -1 for INF bound, and +1 for SUP bound.
    bound_sign = hk.bound(ia)*2-3;
    freeze = hk.freeze(ia);
    
    % Compute the lagrange multipliers oriented wrt the bound direction.
    [l r] = max( -lambda{k} .* bound_sign .* not(freeze) );

    massert( (l==0) | (~hk.freeze(r)),'r is freeze and should not be')
    
    if l>maxl && hk.btype(r)~=Etwin & ~hk.freeze(r)
        maxl=l; cst=[k r];
    end
    
    clear hk;
end

need = (maxl>0);

% --- FREEZE LOOP ----------------------------------------------------
% Instead of keeping all the multipliers of any levels, the algorithm
% "freezes" the constraint corresponding to nonzero multiplier. These
% constraints will never be lexicographic positiv for the following steps of
% the algorithm. A freezed constraint cannot be removed from the active set
% and therefore stays, freezed and active, as an equality constraint until
% the remaining of the active search.
if nargout>=3 && not(need)
    for k=1:kl
        positive                 = find( abs(lambda{k})>THR );
        iPositive                = h(k).active(positive);
        h(k).freeze( iPositive ) = 1;
    end
end

% -----------------------------------------------------------------------
% --- Validity test (not necessary in release mode).
% --- Check the validity of the multipliers ie A'lambda = 0.
if kl<length(h) & norm(active_rows(h(1:kl))'*stacked_cell(lambda))>1e-5
    disp('Error in reconstruction')
    disp(norm(active_rows(h(1:kl))'*stacked_cell(lambda)));
end

