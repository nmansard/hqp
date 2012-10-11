% From the stack matrix Au, build the cell hierarchized matrix A{}.
function A = cut_Au( Au,m );

p=length(m);
A = {};
for k=1:p
    mprec = sum( m(1:k-1) );
    A{k}  = Au( mprec+1:mprec+m(k),: );
end
