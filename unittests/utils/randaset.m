
function [ aset abound  ] = randaset( etype );

n=length(etype);
Etwin = 1;Edouble = 2;Einf = 3;Esup = 4;Enone=0;

for k=1:n
    s=length(etype{k});
    as=floor(rand*(s+1));
    p=randperm(s);
    aset{k} = p(1:as)';

    abound{k}=zeros(0,1);
    for i=1:as
        c=aset{k}(i);
        if etype{k}(c)==Etwin
            abound{k}(i,1) = 1;
        elseif etype{k}(c)==Einf
            abound{k}(i,1) = 1;
        elseif etype{k}(c)==Esup
            abound{k}(i,1) = 2;
        elseif etype{k}(c)==Edouble
            abound{k}(i,1) = round(rand)+1;
        end
    end
end


