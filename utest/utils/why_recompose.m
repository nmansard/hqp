function [ Au ] = why_recompose(h,Y)

% why_recompose recomposes the original stacked matrix from the HCOD
%                 decomposition.
%% Synopsis:
%    Au = why_recompose(h,y)
%% Input:
%    h   the "h" structure built by the HCOD function.
%    Y   the corresponding Y basis.
%% Output:
%    Au  = W*H*Y' the matrix used to generate h.
%
%

p=length(h);
Au=[];
for k=1:p
    hk=h(k);
    iw=hk.iw; im=hk.im; il=hk.im(hk.n+1:hk.m);
    jm=1:hk.rp;       % Indices of the columns corresponding to M and N.
    jl=hk.rp+1:hk.ra; % Indices of the columns corresponding to L.
    
    Au = [Au; ...
              hk.W(iw,im)*hk.H(im,jm)*Y(:,jm)' ...
          +   hk.W(iw,il)*hk.H(il,jl)*Y(:,jl)' ...
          ];
end
