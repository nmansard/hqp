function blabla (a);
str='';
for i=1:length(a)
    stri= num2str(a{i});
    str=strcat(str,stri(1:end));
end
disp(str);
