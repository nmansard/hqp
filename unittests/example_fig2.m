% Compute the optimum of a 2D problem and draws a representation of the
% active-search path with respect to the constraints.
% The problem is the one drawn in the paper in Fig. 2. 
% To comply with the example of the paper, the active-search algo should
% start from a given initial point, which is not possible with the hiearchic
% solver. To produce a similar behavior, the problem is shift, so that x0 is
% at the origin. And a third task is added as "x = -x0" to attract the
% solver to the origin of the non-shifted problem.
% The original problem is first plotted. Then, it is shifted and given to
% the solver. The path of the solver is then plotted.

addpath('unittests/utils');
addpath('utils');
constants;

% Definition of the problem: two levels and an initial point.
clear A b x0;
A{1} = [  0.1 -1 ;  1 -1  ];   b{1} = [ -0.55 ;   1.5   ];
A{2} = [ -1   -1 ; -1  0  ];   b{2} = [ -2    ;  -2.5   ];
x0   = [ -0.5 ; 1.5 ];

% --- CONSTRAINT PLOT ------------------------------------------------
% --- CONSTRAINT PLOT ------------------------------------------------
% --- CONSTRAINT PLOT ------------------------------------------------

% Plot the non shifted problem.

clf;  hold on
color = [     0.75 0.75 0;    0 0 0.78;    1 0  0;   0 0 0  ];

xpoints = -10:1:10;
plot(xpoints,0*xpoints,'color',[.9 .9 .9])
plot(0*xpoints,xpoints,'color',[.9 .9 .9])

for k=1:length(A)
    ck=color(k,:); ck2 = ck/norm(ck)*0.2+[ .8 .8 .8];
    for i=1:size(A{k},1)
        Ai=A{k}(i,:); d=b{k}(i);
        a1=Ai(1); a2=Ai(2);
        if abs(a1)>abs(a2)
            hplot(k) =  plot( (d-a2*xpoints)/a1,xpoints, 'color', ck );
            plot( (d+1e-1-a2*xpoints)/a1,xpoints, 'color', ck2, 'LineWidth',3 );
        else
            hplot(k) = plot( xpoints,(d-a1*xpoints)/a2, 'color', ck);
            plot( xpoints,(d+1e-1-a1*xpoints)/a2, 'color', ck2, 'LineWidth',3 );
        end
    end
end

axis([ -1.5 4  -1 2.5 ]);

% --- SOLVER EXECUTION -----------------------------------------------
% --- SOLVER EXECUTION -----------------------------------------------
% --- SOLVER EXECUTION -----------------------------------------------

% Problem resolution.

% The problem is first shifted to have x0 at the origin.
b{1} = b{1}-A{1}*x0;
b{2} = b{2}-A{2}*x0;
% A third task is added as "x = -x0" to attract the solver to the origin of the
% non-shifted problem.
A{3} = [ 1 0 ; 0 1 ]; b{3} = -x0;

% Write the problem under a form that is appropriate for the solver.
for k=1:length(b) b{k} = [b{k} b{k} ]; end
btype{1} = [ Esup Esup ];
btype{2} = [ Esup Esup ];
btype{3} = [ Etwin Etwin ];
[aset abound] = dummybound(zeros(3,1));

% Problem resolution.
[primal dual h Y xtrack ]=active_search_verbose(A,b,btype,aset,abound);

% Check the proper behavior of the solver.
tracking = [
   0.00000   0.32258   0.44554   1.81818   2.35294   2.77778   2.77778   3.00000
   0.00000  -0.96774  -0.95545  -0.81818  -0.76471  -0.72222  -0.72222  -0.50000
];
massert( norm([xtrack{1,:}]-tracking)<1e-3, ...
         'The result is not what was expected.');

% --- EXEC PLOT ------------------------------------------------------
% --- EXEC PLOT ------------------------------------------------------
% --- EXEC PLOT ------------------------------------------------------

% Plot the algorithm path. 

% The path is first shifted back to fit the original problem.
xi=[xtrack{1,:}](1,:)+x0(1);
yi=[xtrack{1,:}](2,:)+x0(2);

% The initial point is plotted as a circle.
hplot(end+1) = plot(x0(1),x0(2),'o','color',[.6 .6 .6],'MarkerSize',20,'LineWidth',3);
% The path is plotted in red with cross +.
hplot(end+1) = plot(xi,yi,'r-+','LineWidth',2,'MarkerSize',12);
% The final point is plotted as a red star *.
hplot(end+1) = plot(xi(end),yi(end),'r-*','LineWidth',2,'MarkerSize',20);

title('Constraint position and active-search path')
%legend(hplot,'Level 1','Level 2','Initial point','A-search path','Final point');
