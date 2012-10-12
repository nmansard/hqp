function [ h ] = freeze(lambda,h,THR);

% freeze stores a Boolean in "h" to prevent the deactivation of the constraints
% whose multipliers are nonzero. 
%% Synopsis:
%    h = freeze(lambda,h)
%    h = freeze(lambda,h,THR)
%% INPUT:
%    lambda: the multipliers for level kl of h.
%    h       "h" structure storing all the HQP data.
%    THR     threshold to test the positivity.
%% Output:
%    h       The function has a side effect: it modifies the "freeze" field of 
%               the "h" structure. To account for the side effect, h is returned.
%
% Instead of keeping all the multipliers of any levels, the algorithm "freezes"
% the constraint corresponding to nonzero multiplier. These constraints will
% never be lexicographic positiv for the following steps of the algorithm. A
% freezed constraint cannot be removed from the active set and therefore stays,
% freezed and active, as an equality constraint until the remaining of the
% active search.
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- cf. COPYING.LESSER
%
% --- DEFAULT ARGUMENTS --------------------------------------------------------
if nargin==2
    THR=1e-8;
end
% ---------------------------------------------------------------------

kl = length(lambda);

for k=1:kl
    positive                 = find( abs(lambda{k})>THR );
    iPositive                = h(k).active(positive);
    h(k).freeze( iPositive ) = 1;
end

