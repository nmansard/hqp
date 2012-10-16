addpath('unittests/utils');
addpath('utils');
constants;

% --- CONSTRAINT PLOT ------------------------------------------------
% --- CONSTRAINT PLOT ------------------------------------------------
% --- CONSTRAINT PLOT ------------------------------------------------
clf;  hold on
color = [     0.75 0.75 0;    0 0 0.78;    1 0  0;   0 0 0  ];

A={}; b={};

A{1} = [ -1 -0.1 ;  -1 -1];  b{1} = [1 ; 3.5];
A{2} = [  -1 1 ;  0 1   ];  b{2} = [-1; -3   ];
central = [ -1.5 ; -0.5 ];

for k=1:length(A)
    ck=color(k,:); ck2 = ck/norm(ck)*0.2+[ .8 .8 .8];
    for i=1:size(A{k},1)
        Ai=A{k}(i,:); d=b{k}(i);
        a1=Ai(1); a2=Ai(2);
        x=-100:1:100;
        if abs(a1)>abs(a2)
            plot( (d-a2*x)/a1,x, 'color', ck );
            plot( (d+1e-1-a2*x)/a1,x, 'color', ck2, 'LineWidth',3 );
        else
            plot( x,(d-a1*x)/a2, 'color', ck);
            plot( x,(d+1e-1-a1*x)/a2, 'color', ck2, 'LineWidth',3 );
        end
       
    end
end

axis([ -3 3 -5 1]);

plot(-100:100,0*(-100:100),'k')
plot(0*(-100:100),-100:100,'k')
plot(central(1),central(2),'+','color',[.6 .6 .6],'MarkerSize',20,'LineWidth',5);

% --- SOLVER EXECUTION -----------------------------------------------
% --- SOLVER EXECUTION -----------------------------------------------
% --- SOLVER EXECUTION -----------------------------------------------
%A{2} = A{2}([2 1],:);
A{3} = [ 1 0 ; 0 1 ]; b{3} = central;
for k=1:length(b) b{k} = [b{k} b{k} ]; end

btype{1} = [ Esup Esup ];
btype{2} = [ Esup Esup ];
btype{3} = [ Etwin Etwin ];
[aset abound] = dummybound(zeros(3,1));

[primal dual h Y xtrack ]=active_search_verbose(A,b,btype,aset,abound);

% --- EXEC PLOT ------------------------------------------------------
% --- EXEC PLOT ------------------------------------------------------
% --- EXEC PLOT ------------------------------------------------------
plot([xtrack{1,:}](1,:),[xtrack{1,:}](2,:),'r-+','LineWidth',2,'MarkerSize',12)
plot([xtrack{1,:}](1,end),[xtrack{1,:}](2,end),'r-*','LineWidth',2,'MarkerSize',20)

tracking = [
   0.00000  -0.96774  -0.95545  -0.81818  -0.76471  -0.72222  -0.72222  -0.50000
   0.00000  -0.32258  -0.44554  -1.81818  -2.35294  -2.77778  -2.77778  -3.00000
];
massert( norm([xtrack{1,:}]-tracking)<1e-3, ...
         'The result is not what was expected.');
