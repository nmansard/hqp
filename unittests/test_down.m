addpath('unittests/utils');

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[3 4 3 4],[1 3 3 3]);
[active,bound] = dummybound(mref);
[h,Y] = hcod(A,b,btype,active,bound);
A{1}(1,:) = Y(end,:);
[h,Y] = hcod(A,b,btype,active,bound);

% ------------------------------------------------------------------------------
% Remove a row from a rank deficient level (no regularization no promotion).
testname = 'DOWN1-RD: ';
[hd,Yd] = down(1,2,h,Y);
if norm(why_recompose(hd,Yd)-active_rows(hd))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
end
hdr = [hd.r];
hr = [h.r];
if any(hdr(1)~=hr(1))
    blabla({testname,'Error in the rank of the modified level.'});
end
if any(hdr(2:end)~=hr(2:end))
    blabla({testname,'Error in the rank of the other levels.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Remove a full-rank row, leading to a regularization but no promotion.
testname = 'DOWN1-FR: ';
[hd,Yd] = down(1,1,h,Y);
if norm(why_recompose(hd,Yd)-active_rows(hd))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
end
hdr = [hd.r];
hr = [h.r];
if hdr(1)~=hr(1)-1
    blabla({testname,'Error in the rank of the modified level.'});
end
if any(hdr(2:end)~=hr(2:end))
    blabla({testname,'Error in the rank of the other levels.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Remove a full-rank row, leading to regularization and promotion.
testname = 'DOWN1-FRPROM: ';
[hd,Yd] = down(3,1,h,Y);
if norm(why_recompose(hd,Yd)-active_rows(hd))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
end
hdr = [hd.r];
hr = [h.r];
if hdr(3)~=hr(3)-1
    blabla({testname,'Error in the rank of the modified level.'});
end
if any(hdr(1:2)~=hr(1:2))
    blabla({testname,'Error in the rank of the other levels.'});
end
if hdr(4)~=hr(4)+1
    blabla({testname,'Error in the rank of the promoted levels.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Remove all the rows of a level.
testname = 'DOWN1-ALL: ';
[hd,Yd] = down(1,3,h,Y);
[hd,Yd] = down(1,2,hd,Yd);
[hd,Yd] = down(1,1,hd,Yd);
if norm(why_recompose(hd,Yd)-active_rows(hd))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
end
hdr = [hd.r];
hr = [h.r];
if hdr(1)~=0
    blabla({testname,'Error in the rank of the modified level.'});
end
if any(hdr(2:4)~=hr(2:4)+[1 0 0])
    blabla({testname,'Error in the rank of the other levels.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Remove all the rows of a level.

testname = 'DOWN2-DEF: ';
[ A b btype nh p mref rref Au bu ] = randstack(12,4,[2 7 1 2],[2 4 1 1]);
[active,bound] = dummybound({ 1:2, 1:6, 1, 1:2 });
[h,Y] = hcod(A,b,btype,active,bound);

[hd,Yd] = down(2,5,h,Y);
if norm(why_recompose(hd,Yd)-active_rows(hd))> 1e-6
    blabla({testname,'Error in the reconstruction.'});
end
hdr = [hd.r];
if any(hdr~=[2 4 1 1])
    blabla({testname,'Error in the rank.'});
end

%disp('Next'); return
