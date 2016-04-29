function [ temp ] = twoscomp( val1)%,val2 )

%v = typecast([val1, val2], 'uint16');
v=val1;
temp=single(0);

if bitand(v,1)==1
    temp=double(temp+1);
end
if bitand(v,2)==2
    temp=double(temp+2);
end
if bitand(v,4)==4
    temp=double(temp+4);
end
if bitand(v,8)==8
    temp=double(temp+8);
end
if bitand(v,16)==16
    temp=double(temp+16);
end
if bitand(v,32)==32
    temp=double(temp+32);
end
if bitand(v,64)==64
    temp=double(temp+64);
end
if bitand(v,128)==128
    temp=double(temp+128);
end
if bitand(v,256)==256
    temp=double(temp+256);
end
if bitand(v,512)==512
    temp=double(temp+512);
end
if bitand(v,1024)==1024
    temp=double(temp+1024);
end
if bitand(v,2048)==2048
    temp=double(temp+2048);
end
if bitand(v,4096)==4096
    temp=double(temp+4096);
end
if bitand(v,8192)==8192
    temp=double(temp+8192);
end
if bitand(v,16384)==16384
    temp=double(temp+16384);
end
if bitand(v,32768)==32768
    temp=double(temp-32768);
end