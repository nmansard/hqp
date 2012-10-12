addpath('utest/utils');

% Compute the solution of the eHQP using the iterative equations of
% Siciliano91.
function x = siciliano( Aact,bact,m,EPS );
   
    Aact = cut_Au(Aact,m);
    bact = cut_Au(bact,m);
    nh   = columns(Aact{1});
    p    = length(Aact);
    
    x    = zeros(nh,1);
    P    = eye(nh);

    for k=1:p
        AP = (Aact{k}*P);
        AP_inv = pinv(AP,EPS);
        x = x + AP_inv*( bact{k} - Aact{k}*x );
        P = P - AP_inv*AP;
    end
end


% ------------------------------------------------------------------------------
% Compare for a 1-level full-rank fully-active stack.
% The HQP is full-row rank, the h-inverse is the pseudoinverse.
testname = 'SOLVE-1L: ';

[ A b btype nh p mref rref Au bu ] = randstack(6,1,3,3);
[active,bound] = dummybound( { 1:3 } );
[h,Y] = hcod(A,b,btype,active,bound);
[y x] = ehqp_primal(h,Y);

[Aact,bact] = active_rows(h) ;
x_pinv = pinv(Aact)*bact;

if norm(x-x_pinv)> 1e-6
    blabla({testname,'Error in the norm.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Compare for a 2-level full-rank fully-active stack
testname = 'SOLVE-2L: ';

[ A b btype nh p mref rref Au bu ] = randstack(6,2,[ 2 3 ],[ 2 3 ]);
[active,bound] = dummybound( { 1:2 1:3 } );
[h,Y] = hcod(A,b,btype,active,bound);
[y x] = ehqp_primal(h,Y);

[Aact,bact] = active_rows(h) ;
x_pinv = pinv(Aact)*bact;

if norm(x-x_pinv)> 1e-6
    blabla({testname,'Error in the norm.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Compare for a partially activated multi level stack
testname = 'SOLVE 4L-FR: ';

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[3 4 3 4],[2 3 3 3]);
[active,bound] = dummybound( { 1, 1:3, 1:2, 3:4 });
[h,Y] = hcod(A,b,btype,active,bound);
[y x] = ehqp_primal(h,Y);

[Aact,bact] = active_rows(h) ;
x_pinv = pinv(Aact)*bact;

if norm(x-x_pinv)> 1e-6
    blabla({testname,'Error in the norm.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Compare for a 2-level rank-def fully-active stack
testname = 'SOLVE-3L-RD: ';

[ A b btype nh p mref rref Au bu ] = randstack(6,3,[ 4 4 4],[ 2 3 1 ]);
[active,bound] = dummybound( { 1:4, 1:4, 1:4, 1:4 });
[h,Y] = hcod(A,b,btype,active,bound);
[y x] = ehqp_primal(h,Y);

[Aact,bact] = active_rows(h) ;
x_sici = siciliano(Aact,bact,[h.m],1e-5);
if norm(x-x_sici)> 1e-5
    blabla({testname,'Error in the norm.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------
% Compare for a partially activated rank-def multi level stack
testname = 'SOLVE 4L-RD: ';

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[5 4 6 4],[2 2 3 3]);
[active,bound] = dummybound({ [ 1 3 5 ], 1:3, [ 1 2 3 5 ],  1:2 });
[h,Y] = hcod(A,b,btype,active,bound);
[y x] = ehqp_primal(h,Y);

[Aact,bact] = active_rows(h) ;
x_sici = siciliano(Aact,bact,[h.m],1e-5);
if norm(x-x_sici)> 1e-5
    blabla({testname,'Error in the norm.'});
end

%disp('Next'); return
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% ---HT SOLVE ------------------------------------------------------------------
% ------------------------------------------------------------------------------

% Compare for a partially activated rank-def multi level stack
testname = 'TSOLVE 4L-RD: ';

[ A b btype nh p mref rref Au bu ] = randstack(12,4,[5 4 6 4],[2 2 3 3]);
[active,bound] = dummybound({ [ 1 3 5 ], 1:3, [ 1 2 3 5 ],  1:4 });
[h,Y] = hcod(A,b,btype,active,bound);
[y x] = ehqp_primal(h,Y);

k=4;
bk=_b(h,k);
Ak=_A(h,k);

rhok = Y'*Ak'*(Ak*x-bk);
lambda = ehqp_dual(k,y,h,Y);

if norm(active_rows(h(1:k))' * stacked_cell(lambda)) > 1e-5
    blabla({testname,'Error in the norm of C^T lambda.'});
end

%disp('Next'); return

% ------------------------------------------------------------
testname = 'TSOLVE 4L-RD LAST: ';
mult = ehqp_dual(5,y,h,Y);

if norm(active_rows(h)' * stacked_cell(mult) + x) > 1e-5
    blabla({testname,'Error in the norm of last stage C^T lambda.'});
end

%disp('Next'); return
