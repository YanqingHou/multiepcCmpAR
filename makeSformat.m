function str=makeSformat(n)
if n<=0
    error('wrong number of s!');
end
str='\n';
for i=1:n
    str=strcat('%12.5f ',str);
end

