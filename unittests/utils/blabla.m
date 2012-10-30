function blabla (a);
str='';
for i=1:length(a)
    stri= num2str(a{i});
    if strcmp(stri,'\n')
        stri=char(10);
    end
    str=strcat(str,stri(1:end));
end
disp(str);
