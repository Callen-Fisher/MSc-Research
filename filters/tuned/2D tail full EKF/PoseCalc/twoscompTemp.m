function [ temp ] = twoscompTemp( val1 )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

v= val1;

temp=single(0);
if((bitand(v,128))==128)
    temp=single(-128);
end
if((bitand(v,1))==1)
    temp=single(temp+1);
end
if((bitand(v,2))==2)
    temp=single(temp+2);
end
if((bitand(v,4))==4)
    temp=single(temp+4);
end
if((bitand(v,8))==8)
    temp=single(temp+8);
end
if((bitand(v,16))==16)
    temp=single(temp+16);
end
if((bitand(v,32))==32)
    temp=single(temp+32);
end
if((bitand(v,64))==64)
    temp=single(temp+64);
end

end
