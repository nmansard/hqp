% Display nicely the matrix A.
function pm(A,name)

if nargin==2
    disp(sprintf('%s = ',name));end
[n m]=size(A);

if( sumsq(sumsq(A))==Inf )
  f=1;
else
  f=(round(log(norm(A))/log(10)/3)*3);
end
if(f!=0) disp(sprintf('1e%d x',f)); end
f=10^f;
for i=1:n
  s='';
  for j=1:m
    a=A(i,j);
    if(abs(a)<1e-6*f) s=sprintf('%s   0       ',s);
    else if(abs(a)<1e-3*f) s=sprintf('%s   0...    ',s);
    else
      if(a<0) sl=sprintf('  -'); else sl=sprintf('   '); end;
      sl=sprintf('%s%.5f',sl,abs(a/f));
      while length(sl)<11 sl=sprintf('%s ',sl); end
      s=sprintf('%s%s',s,sl);
    end end
  end
  disp(s);
end
