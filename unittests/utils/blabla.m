function blabla (a);
str='';
for i=1:length(a)
    stri= disp(a{i});
    str=strcat(str,stri(1:end-1));
end
disp(str);
