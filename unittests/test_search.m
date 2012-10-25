addpath('unittests/utils');

seed=round(rand*10000);
%seed = 283;8394
%seed = 5123;
%seed = 8394;

[ nh p m r ] = randstackref(seed);
blabla({'seed = ',seed});
blabla({'size = ',nh,' x',p});
[ A b btype nh p mref rref Au bu ] = randstack(nh,p,m,r);
%[aset abound] = randaset(btype);
[aset abound] = dummybound(zeros(p,1));

[primal dual h Y ]=active_search_verbose(A,b,btype,aset,abound);

blabla({'Primal optimum = \n',primal'});
blabla({'\nOptimum active set =\n',dispaset(h)});

%hdist(primal,h)

    

