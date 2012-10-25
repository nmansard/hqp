addpath('unittests/utils');

% Unitary test to check the validity of the hcod.m function. Decompose a
% random matrix Au and check that the recomposed decomposition is equal to
% the original matrix.

[ A b btype nh p mref rref Au bu ] = randstack(12);
for k=1:p
    active{k}=(1:size(A{k},1))';
    bounds{k} = ones(length(active{k}),1);
end

[h,Y] = hcod(A,b,btype,active,bounds);

for k=1:p
    if h(k).r ~= rref(k)
        blabla({'Error in the rank of the decomposition at level ',k});
    end
end
   
Acheck=cut_Au(why_recompose(h,Y),mref);
for k=1:p
    if norm(Acheck{k}-A{k})>1e-6
        blabla({'Error in the reconstruction at level ',k});
    end
end


