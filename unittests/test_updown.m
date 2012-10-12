addpath('utest/utils');

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[3 4 3 4],[2 3 3 3]);
[active,bound] = dummybound( { 1, 1:3, 1:2, 3:4 } );
[h,Y] = hcod(A,b,btype,active,bound);

% ------------------------------------------------------------------------------
% Introduce a full rank Ju at the last stage.
testname = 'UP4-FR +4.1-4.1: ';
[hu,Yu] = up(4,1,'right',h,Y);
[hd,Yd] = down(4,hu(4).m,hu,Yu);

if norm(why_recompose(hu,Yu)-active_rows(hu))> 1e-6
    blabla({testname,'Error in the u reconstruction.'});
end
if norm(why_recompose(hd,Yd)-active_rows(hd))> 1e-6
    blabla({testname,'Error in the d reconstruction.'});
end
if norm(why_recompose(hd,Yd)-active_rows(h))> 1e-6
    blabla({testname,'Error in the +u-d reconstruction.'});
end
if any([hu.r]!=([h.r]+[0 0 0 1]))
    blabla({testname,'Error in the rank.'});
end
if any([hd.r]!=([h.r]))
    blabla({testname,'Error in the rank of th +u-d.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Introduce a full rank Aup at the last stage.
testname = 'UP4-FR -1.1+1.1: ';
k=1; r=1; 
[hd,Yd] = down(k,r,h,Y);
[hu,Yu] = up(k,h(k).active(r),'left',hd,Yd);

if norm(why_recompose(hu,Yu)-active_rows(hu))> 1e-6
    blabla({testname,'Error in the u reconstruction.'});
end
if norm(why_recompose(hd,Yd)-active_rows(hd))> 1e-6
    blabla({testname,'Error in the d reconstruction.'});
end
c0=sum([h(1:k-1).m])+r; c1=sum([h(1:k).m]);
S=eye(sum([h.m]));
S(c0:c1,:) = S([c0+1:c1 c0],:);
if norm(why_recompose(hu,Yu)-S*active_rows(h))> 1e-6
    blabla({testname,'Error in the +u-d reconstruction.'});
end
if any([hu.r]!=([h.r]))
    blabla({testname,'Error in the rank of th +u-d.'});
end

%disp('Next'); return


