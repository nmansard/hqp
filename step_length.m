function [ need taumax cst] = step_length(x0,x1,h,Y,THR);

% step_length  checks if the current solution x1 respects the unactive
%                  constraints.  In that case, return [true 1 []].
%                  Otherwise, return false, the step length and the
%                  reference on the first violated constraint. This function
%                  corresponds roughly to Alg. 4 of the paper.
%% Synopsis:
%    [ need taumax cst ] = step_length(x0,x1,h,Y)
%    [ need taumax cst ] = step_length(x0,x1,h,Y,THR)
%% Input:
%    x0      current value of the primal optimum.
%    x1      next value of the primal optimum.
%    h       "h" structure storing all the HQP data.
%    Y       right basis of the HCOD.
%    THR     is the threshold used to test the collision.
%% Output:
%    need    is true iff at least one constraint needs to be activated at
%               the next step of the active search.
%    taumax  is the step length.
%    cst     is the constraint reference if one constraint needs to be
%               activated. In that case, cst = [ k c bound ], the three
%               arguments to be given to the "up" function.
%
% The function does not follows exactly Alg. 4 because of the double
% bounds. It uses the subfunctions check_bound to factorize the code of
% Alg. 4#6 and 4#10. Moreover, it uses a shortcut by setting a step of
% length tau=1-EPS if a bound is not satisfied for x_0.
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- cf. COPYING.LESSER
%

% --- DEFAULT ARGUMENTS --------------------------------------------------------
nin = nargin
if nin==4
    % Default argument for the HCOD threshold.
    THR=1e-8;
    nin=nin+1;
end
% ---------------------------------------------------------------------

p=length(h);
nh=size(Y,2);
constants;

taumax=1; cst=[];
for k=1:p
    massert( not(exist('hk')) ,'Error, hk should have been cleared.');
    hk=h(k);
    mmax=hk.mmax;

    for r=setdiff(1:mmax,hk.active);
        b   = hk.b(r,:);
        typ = hk.btype(r);
        Ax1 = hk.A(r,:)*x1;
        
        [ typ1 b1 ] = check_bound(Ax1,b,typ,THR);  % Alg 4#10

        if typ1~=Enone
            Ax0  = hk.A(r,:)*x0;
            typ0 = check_bound(Ax0,b,typ,THR);     % Alg 4#6
            if typ0==Enone
                tau=(b1-Ax0)/(Ax1-Ax0);            % Alg 4#7
            else
                tau=1-THR;                         % Alg 4#9
            end

            if tau<taumax                          % Alg 4#15
                taumax=tau;
                if typ1==Einf cst=[k r 1];
                else          cst=[k r 2];
                end
            end
        end
    end
    clear hk;
end

need=taumax<1;                                    % Alg 4#17 and 4#11

end

% --- SUBFUNCTIONS --------------------------------------------------

function [ violtype violval ] = check_bound(val,bound,type,THR);

% check_bound   checks if the input value "val" is inside the bounds specified 
%                 by bound,type.
%% Synopsis:
%   violtype             = check_bound(val,bound,type)
%   [ violtype violval ] = check_bound(val,bound,type);
%   [ violtype violval ] = check_bound(val,bound,type,THR);
%% Input:
%    val       is the value to test against bound.
%    bound     is a 1x2 vector of floats.
%    type      specifies the bound type (=,<=,>=,<=<=).
%    THR       is the threshold used to test the collision.
%% Output:
%    violtype  if no bound is violated, return Enone, otherwise returns Einf or
%                Esup depending on the violation.
%    violval   is the value of "bound" that is met.
%
%

    constants;

    violtype=Enone; violval=0;
    if (type==Einf | type==Edouble) & (val<bound(1)-THR)
        violtype=Einf; violval=bound(1);
    end
    if (type==Esup | type==Edouble) & (val>bound(2)+THR)
        violtype=Esup; violval=bound(2);
    end

end


