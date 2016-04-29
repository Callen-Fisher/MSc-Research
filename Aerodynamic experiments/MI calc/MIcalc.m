clc
clear




disBetweenIVdisc=[24,30,24,34,34,42,49,52,50,50,48,44,40,37,32,30,25,22,22,16,17]*1E-3;%length of bone 
disCranialEndIV=[2.8,6.3,8.8,12,15.8,20,25,32,35.2,40,45,49.2,53,56.7,60,62.9,65.6,68,70.4,72.3,74]*1E-2;
individualMass=[50.8,61.2,56,53.3,55.2,48.4,45.9,32.2,28.7,20.4,15.1,11.8,8,5.6,4.2,2.9,2.2,1.6,1,0.7,0.3]*1E-3;%+0.162/21;
diameterOfBone=[19,11,17,13,14,14,13,12,13,13,12,10,9,9,8,7,6,6,6,6,5]*1E-3;
diameterOfBoneAndMuscle=[70,90,71,65,56,46,38,27,25,23,24,18,17,16,12,12,10,9,7,6,5]*1E-3;
mean(diameterOfBoneAndMuscle)
sum(disBetweenIVdisc)
%plot(cumsum(disBetweenIVdisc),disCranialEndIV)


Ix=zeros(1,length(disBetweenIVdisc));

%individual moment of inertia of the bones
for i=1:1:length(disBetweenIVdisc)
    Ix(i)=(1/12)*individualMass(i)*(3*(diameterOfBoneAndMuscle(i)*0.5)^2+disBetweenIVdisc(i)^2);
end

%net moment of inertia
netMoment=0;
for i=1:1:length(disBetweenIVdisc)
    netMoment=netMoment+Ix(i)+individualMass(i)*(disCranialEndIV(i)-0.5*disBetweenIVdisc(i))^2;% double check the H^2 term
    temp=disCranialEndIV(i)-0.5*disBetweenIVdisc(i);
end
netMoment=netMoment+1/12*0.162*(3*((0.03/2)^2+(0.023/2)^2)+0.74^2)+0.162*(0.74/2)^2;%adds the mass of the tail
netMoment

%COM
netTorque=0;
for i=1:1:length(disBetweenIVdisc)
    netTorque=netTorque+individualMass(i)*(disCranialEndIV(i)-0.5*disBetweenIVdisc(i));%double check the H term
end
COM=netTorque/sum(individualMass)% the sum of the mass is only 0.5KG----take into account the fur? 
sum(individualMass);


comPercentage=COM/sum(disBetweenIVdisc)
