function [ lambdak ]= ehqp_dual(kl,y,h,Y);

%  ehqp_dual   computes the dual optimum (Lagrange multipliers for stage "kl",
%                 being given the primal optimum expressed in the Y basis. 
%                 This function corresponds to Alg. 2 of the paper. The w_k^* 
%                 is also computed.
%% Synopsis:
%    lambda   = ehqp_dual(kl,y,h,Y)
%% Input:
%    kl      id of the level for which the multipliers are computed.
%    y       is the primal optimum in the Y basis computed by ehqp_primal.
%    h       "h" structure storing all the HQP data.
%    Y       right basis of the HCOD.
%% Output:
%    lambdak multipliers corresponding to the level "kl". The results is
%              partitioned for each levels i=1:kl into a cell. w_kl^*
%              corresponds to lambdak{kl}.
%
% This function corresponds to Alg. 2 of the paper "Hiearchical Quadratic
% Programming" (Part I).
%
% Copyright Nicolas Mansard -- LAAS/CNRS
%    -- and Adrien Escande -- JRL/CNRS
%    -- cf. COPYING.LESSER
%

p=length(h);
nh=size(Y,2);
lambdak = {};

% Octave needs to know that y is a col vector.
if length(y)==1 y(2,1)=0; end

% ------------------------------------------------------------------------------
% 1. Compute rho: rho = Y'*J'*( e-J*x )
if kl>p
    % 1. [Level p+1]: Try to minimize the norm of the solution: this case is explored
    % in Sec. I-4.5.
    rhobar = -y;
else

    % 1. [Level k<=p]: Compute the slack w in the W basis: 
    %         ===>  W'*w = W'*A*x - b = M0*ybar - W0'b
    massert( not(exist('hk')) ,'Error, hk should have been cleared.');
    hk=h(kl);
    iw=hk.iw; im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;

    im0      = im(1:n);
    im1      = im(n+1:end);
    M0       = hk.H(im0,1:rp);
    W0       = hk.W(iw,im0);
    b        = hk.b( hk.activeb );
    if size(b,2)==0 b=b'; end % Octave patch
    
    rho         =  M0*y(1:rp) - W0'*b;  % Alg 2#4
    lambdak{kl} =  W0*rho;              % Alg 2#5
    rhobar      = -M0'*rho;             % Alg 2#6
    
    clear hk;
    clear iw im r n ra rp m im0 im1 M0 W0 b;
end
 
% ------------------------------------------------------------------------------
% 2.a Compute for the first kl-1 levels the multipliers: lambda = W' [ML]^H  rho

for k=kl-1:-1:1
    massert( not(exist('hk')) ,'Error, hk should have been cleared.');
    hk=h(k);
    iw =hk.iw;   im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;
    im0=im(1:n); im1=im(n+1:end);

    if r>0
        L          = hk.H(im1,rp+1:ra);
        M1         = hk.H(im1,1:rp);
        W1         = hk.W(iw,im1);

        rho        = rhobar(rp+1:ra);  % Alg 2#8
        rhobar     = rhobar(1:rp,1);
    
        rho        = L'\rho;           % Alg 2#9
        lambdak{k} = W1*rho;           % Alg 2#10
        if k>1
            rhobar = rhobar - M1'*rho; % Alg 2#11
        end
    else
        lambdak{k} = zeros(m,1);
    end
    
    clear hk;
    clear iw im r n ra rp m im0 im1 M1 W1 ;
end



