addpath('unittests/utils');

seed=round(rand*10000);

[ nh p m r ] = randstackref(seed);
blabla({'seed = ',seed});
blabla({'size = ',nh,' x',p});
[ A b btype nh p mref rref Au bu ] = randstack(nh,p,m,r);
[aset abound] = dummybound(zeros(p,1));

[primal dual h Y ]=active_search_verbose(A,b,btype,aset,abound);

blabla({'\n','Primal optimum = ','\n',primal'});
blabla({'\n','Optimum active set =','\n',dispaset(h)});

% Display a more complex diagnostic about the solution.
%hdist(primal,h)

    

