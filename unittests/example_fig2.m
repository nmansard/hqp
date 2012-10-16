clf;  hold on
color = [     0.75 0.75 0;    0 0 0.4;    1 0  0;   0 0 0  ];

A={}; b={};

A{1} = [ 1 0.1 ;  1 1];  b{1} = [-1 ; -3.5];
A{2} = [   0 1 ; 1 -1   ];  b{2} = [ -3 ; 1 ];
A{2} = [   0 1 ; 1 -1 ; 1 -2  ];  b{2} = [ -3 ; 1 ; 4];
%A{3} = [ 1 -1  ; -0.1 1];  b{3} = [0; -1];
A{4} = [ 1 0 ;  0 1 ];  b{4} = [0 ; 0];

for k=1:length(A)
    for i=1:size(A{k},1)
        Ai=A{k}(i,:); d=b{k}(i);
        a1=Ai(1); a2=Ai(2);
        x=-100:10:100;
        %if abs(a1)>abs(a2)
        %    plot( -a2/a1*x+d,x, 'color', color(k,:) );
        %else
        %    plot( x,-a1/a2*x+d, 'color', color(k,:) );
        %end
        if abs(a1)>0
            plot( (d-a2*x)/a1,x, 'color', color(k,:) );
        end
        if abs(a2)>0
             plot( x,(d-a1*x)/a2, 'color', color(k,:) );
        end
       
    end
end

axis([ -3 3 -5 1])

plot(-1.111,-1.111,'+','color',[.6 .6 .6],'MarkerSize',20,'LineWidth',5);

x = [ 0 -0.90909 -0.89109 -0.76190 -0.81818 -0.72222 -0.50000];
y = [ 0 -0.90909 -1.08911 -2.38095 -1.81818 -2.77778 -3.00000];

plot( x,y,'r-+','LineWidth',3,'MarkerSize',20)