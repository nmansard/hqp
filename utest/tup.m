addpath('utest/utils');

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[3 4 3 4],[2 3 3 3]);
active = { 1, (1:3)', (1:2)', (3:4)' };
for k=1:p bound{k} = ones(length(active{k}),1); end
[h,Y] = hcod(A,b,btype,active,bound);

% ------------------------------------------------------------------------------
% Introduce a full rank Aup at the last level.
testname = 'UP4-FR: ';
[hup,Yup] = up(4,1,'right',h,Y);
if norm(why_recompose(hup,Yup)-active_rows(hup))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
end
if any([hup.r]!=([h.r]+[0 0 0 1]))
    blabla({testname,'Error in the rank.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Introduce a rank def Ju that does not modify L.
testname = 'UP2-RDM: ';
Adef = h(2).A(4,:);
h(2).A(4,:) = rand(1,3)*h(2).A(1:3,:);
[hu,Yup] = up(2,4,'right',h,Y);
if norm(why_recompose(hu,Yup)-active_rows(hu))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
    pm(why_recompose(hu,Yup)-active_rows(hu));
end
if any([hu.r]!=[h.r])
    blabla({testname,'Error in the rank.'});
end
h(2).A(4,:) = Adef;

%disp('Next'); return
% ------------------------------------------------------------------------------
% Introduce a rank def Ju that modifies L.
testname = 'UP2-RDL: ';
Adef = h(2).A(4,:);
h(2).A(4,:) = rand(1,h(1).m)*h(1).A(h(1).active,:);
[hup,Yupp] = up(2,4,'right',h,Y);
if norm(why_recompose(hup,Yupp)-active_rows(hup))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
    pm(why_recompose(hup,Yupp)-active_rows(hup));
end
if any([hup.r]!=[h.r])
    blabla({testname,'Error in the rank.'});
end
h(2).A(4,:) = Adef;

%disp('Next'); return
% ------------------------------------------------------------------------------
% Introduce a full rank Ju at the first level.
testname = 'UP1-FR: ';
[hup,Yupp] = up(1,2,'right',h,Y);
if norm(why_recompose(hup,Yupp)-active_rows(hup))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
    pm(why_recompose(hup,Yupp)-active_rows(hup));
end
if any([hup.r]!=[h.r]+[1 0 0 0])
    blabla({testname,'Error in the rank.'});
end

%blabla({'Next'}); return
% ------------------------------------------------------------------------------
% Introduce a row deficient wrt the next level.
testname = 'UP1-DEF2: ';
[hup,Yup] = up(2,4,'right',h,Y);
[hup,Yup] = up(1,2,'right',h,Y);
if norm(why_recompose(hup,Yup)-active_rows(hup))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
    pm(why_recompose(hup,Yup)-active_rows(hup));
end
if any([hup.r]!=[h.r]+[1 0 0 0])
    blabla({testname,'Error in the rank.'});
end

%blabla({'Next'}); return
% ------------------------------------------------------------------------------
% Introduce three rows that flush the last level.
clear all
testname = 'UP111-FLUSH: ';

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[3 4 3 6],[3 3 3 3]);
active = { zeros(0,1), (1:3)', (1:3)', (1:6)' };
clear bound; for k=1:p bound{k} = ones(length(active{k}),1); end
[h,Y] = hcod(A,b,btype,active,bound);

[hup,Yup] = up(1,1,'right',h,Y);
[hup,Yup] = up(1,2,'right',hup,Yup);
[hup,Yup] = up(1,3,'right',hup,Yup);
if norm(why_recompose(hup,Yup)-active_rows(hup))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
    pm(why_recompose(hup,Yup)-active_rows(hup));
end
if any([hup.r]!=[h.r]+[3 0 0 -3])
    blabla({testname,'Error in the rank.'});
end

%blabla({'Next'}); return
% ------------------------------------------------------------------------------
% Empty the last level.
testname = 'UP111-EMPTY: ';

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[3 6 3 6],[3 6 3 0]);
active = { zeros(0,1), (1:6)', (1:3)', (1:6)' };
clear bound; for k=1:p bound{k} = ones(length(active{k}),1); end
[h,Y] = hcod(A,b,btype,active,bound);

[hup,Yup] = up(1,1,'right',h,Y);
[hup,Yup] = up(1,2,'right',hup,Yup);
[hup,Yup] = up(1,3,'right',hup,Yup);
if norm(why_recompose(hup,Yup)-active_rows(hup))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
    pm(why_recompose(hup,Yup)-active_rows(hup));
end
if any([hup.r]!=[h.r]+[3 0 0 -3])
    blabla({testname,'Error in the rank.'});
end
if hup(4).r!=0
    blabla({testname,'Error in the rank of the last level.'});
end

%blabla({'Next'}); return
