% Generate a set of non free hierarchical constraint.

if( not(exist('reshot'))) reshot=1; end
if reshot
  clear m n r s J e Ja JP P up;

  % Size of the problem
  m=20;
  % Number of constraints
  n=5;
  % Size of each constraint
  s=[5 5 5 5 8];
  %s=[3 3 3 3];
  % Rank of each constraint
  r=[ 3 5 3 5 3];
  %r=[ 3 3 3 3 ];

  Ja=[]; Pa=eye(m);
  for i=1:n
    Xhi=((-.5+sqrt(rand(s(i),r(i)+size(Ja,1))) )*2);
    J{i}=Xhi*[ ((-.5+rand(r(i),m))*2); Ja ];
    J{i}=J{i}/norm(J{i})*(.75+rand/2);
    e{i}=((-.5+rand(s(i),1))*2);
    %J{i}=J{i}+(rand(s(i),m)-.5)*2*1e-4; % singularity distance

    if(i==3)
      %  J{i}(1,:)=J{i}(4,:);J{i}(2,:)=J{i}(4,:);
      Xhi=((-.5+sqrt(rand(3,1+size(Ja,1))) )*2);
      J{i}=[  Xhi*[ ((-.5+rand(1,m))*2); Ja ]; [ ((-.5+rand(2,m))*2) ] ];
      J{i}=J{i}/norm(J{i})*(.75+rand/2);
    end
    Ja=[Ja;J{i}];
    JP{i}=J{i}*Pa; P{i}=eye(m)-pinv(Ja)*Ja; Pa=P{i};
    %if i==1 upinv{i}=pinv(J{i})*e{i};
    %else upinv{i}=upinv{i-1}+pinv(JP{i})*(e{i}-J{i}*upinv{i-1});end;

    THR=5e-15;
    if i==1 upinv{i}=pinv(J{i},THR)*e{i};
    else upinv{i}=upinv{i-1}+pinv(JP{i},THR)*(e{i}-J{i}*upinv{i-1});end;
  end
end


% --- Test up and down ---
if reshot
  %Jup=rand(1,m)*2-1; k=1;
  %k=3; xhi=rand(1,s(k+1));  Jup=xhi*J{k+1};
  %Jup=J{3}(2,:); k=2;
  k=3; xhi=rand(1,s(k));  Jup=xhi*J{k};
end

%reshot=1;
%reshot=0;

%k=3; ldown=4;
%J{k}(1,:)=J{k}(4,:);
%J{k}(2,:)=J{k}(4,:);
% / --- Test up and down ---
reshot=0;

% --- Reference LQ decomposition
%
% Y: slice of the matrix corresponding to L0 (L*Yp)
% Yp: beginning of the matrix Ya until L0 ([M L]*Yp)
% Z: orthogonal to Yp at iter k, thus not really useful, except Z{n}.
% U,V,W: left basis, W=[V U], W.[N 0 ;M L]=[VN+UM; UL].
% L0: reducted triangle - L: complete right size: L=[ [0;L0] 0 ].
% M: complete left side of the COD: W.[M L].Y'. 
% r: rank of stage k (r(k)=size(L0{k})). s: size of stage k (s=size(J{k},1)).
% ra: used rank at stage k: ra(k)=ra(k-1)+r(k).
clear Y Yp Z U V W L L0 r ra s M;
La=[];
Ya=eye(m);rak=0;
for i=1:n
  Ji=J{i}; si=size(Ji,1);
  JiYa=Ji*Ya; JiY=JiYa(:,1:rak); JiZ=JiYa(:,rak+1:end);
  [Wi Li Qi ri]=cod(JiZ);
  JiY=Wi'*JiY;
  Yqi=eye(m); Yqi(rak+1:end,rak+1:end)=Qi;
  Ya=Ya*Yqi; rak=rak+ri;

  %Yp{i}=Ya(:,1:size(JiY,2));
  %Y{i}=Yi; Z{i}=Zi;
  W{i}=Wi; V{i}=Wi(:,1:si-ri); U{i}=Wi(:,si-ri+1:end);
  L0{i}=Li(1+size(Li,1)-ri:end,1:ri); M{i}=JiY; % L{i}=Li;
  r(i)=ri; if(i==1) ra(i)=ri; else ra(i)=ra(i-1)+ri; end
  s(i)=si;
  %La=[La;JiY Li];
end


if 0%debug
  Yn = [Yp{n} Y{n}]; Zn=Z{n};
  Ya=[Yn Zn];
  Wa=zeros(sum(s));
  for i=1:n
    sprec=sum(s(1:i-1));
    Wa(sprec+1:sprec+s(i),sprec+1:sprec+s(i))=W{i};
  end
end

% check check
% pm(Wa*La*Ya'-Ja)

