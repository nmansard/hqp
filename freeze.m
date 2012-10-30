function [ h ] = freeze(lambda,h,THR);

% freeze stores a Boolean in "h" to prevent the deactivation at the next levels 
% of the constraints whose Lagrange multipliers are nonzero at this level. 
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
% the constraint corresponding to nonzero multipliers. These constraints will
% never be lexicographicaly positive for the following steps of the algorithm. 
% A freezed constraint cannot be removed from the active set and therefore 
% stays, freezed and active, as an equality constraint until the end of the 
% active search.
%
% Copyright Nicolas Mansard -- LAAS/CNRS
%    -- and Adrien Escande -- JRL/CNRS
%    -- cf. COPYING.LESSER
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

