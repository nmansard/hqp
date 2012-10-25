function [ y x ]= ehqp_primal(h,Y);

% ehqp_primal   returns the primal optimum x_p^* and its expression in the Y
%                basis y = Y'*x.  This function follows Alg. 1 of the paper.
%                The slack w^* is not computed here (but in ehqp_dual).
%% Synopsis:
%    y     = ehqp_primal(h,Y)
%    [y x] = ehqp_primal(h,Y)
%% Input:
%    h       "h" structure storing all the HQP data.
%    Y       right basis of the HCOD.
%% Output:
%    y       primal optimum in the Y basis.
%    x       primal optimum in the canonical basis.
%
% Copyright Nicolas Mansard -- LAAS/CNRS -- cf. COPYING.LESSER
%

p=length(h);
nh=size(Y,2);
y = zeros(0,1);

for k=1:p
    massert( not(exist('hk')) ,'Error, hk should have been cleared.');
    hk=h(k);
    iw=hk.iw; im=hk.im; r=hk.r; n=hk.n; ra=hk.ra; rp=hk.rp; m=hk.m;
    if r>0
        im1 = im(n+1:m);
        L   = hk.H(im1,rp+1:ra);
        M1  = hk.H(im1,1:rp);
        W1  = hk.W(iw,im1);
        b   = hk.b( hk.activeb );

        e   = W1'*b - M1*y;        % Alg 1#6
        y(rp+1:ra,1) = L\e;        % Alg 1#8
        %y(rp+1:ra,1) = inv(L)*e;
    end

    h(k)=hk; clear hk;
    clear iw im r n ra rp m im1 W1 M1;
end

if nargout == 2
    x = Y(:,1:h(end).ra)*y;       % Alg 1#11
end

