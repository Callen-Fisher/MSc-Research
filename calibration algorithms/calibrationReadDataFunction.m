function [ gyroS1, magS1, accS1, tempS1, tempADCS1, yawEncoder,yawVelocity] = calibrationReadDataFunction( fileName )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N=1;%the number of sensors that were data logged (gimbal logs 1 and app logs 4)
gimbal=1;
successByte=0;% is there a success byte appended to the packet (tells which sensors were dropped

if(N==4)%if N=4 then using the collar to log so there will be a success byte appended 
    successByte=1;
end

if(gimbal==1)% gimbal always data loggs one sensor
    N=1;
end

fclose('all');                
fid = fopen(fileName);%change for you data file 

rawPackets=[];
%te=[];te=[te; [fread(fid,1,'uint8')]];
while ~feof(fid)%loops through all the data:
  temp=uint8(fread(fid,1,'uint8'));
  if(temp==126)%0x7E
      opcode=fread(fid,1,'uint16');
      if(opcode==257)%0x0101
          tempPacket=[];
          tempRawDataForCRC=[1 1];
          inByte=uint8(fread(fid,1,'uint8'));
          tempRawDataForCRC=[tempRawDataForCRC inByte];
          while(inByte~=126)
              if(inByte==125)%0x7D
                  inByte=uint8(fread(fid,1,'uint8'));
                  tempRawDataForCRC=[tempRawDataForCRC inByte];
                  if(inByte==94)%0x5E
                      tempPacket=[tempPacket 126];%0x7E
                  elseif(inByte==93)%0x5D
                      tempPacket=[tempPacket 125];%0x7D
                  else
                      %return error
                  end
              else
                  tempPacket=[tempPacket inByte];
              end
              inByte=uint8(fread(fid,1,'uint8'));
              tempRawDataForCRC=[tempRawDataForCRC inByte];
          end
          %use tempRawDataForCRC to calc the CRC
          %check CRC here. If pass add the packet, if not then remove
          %packet and raise a flag to show packet was lost
          %CRC=fread(fid,1,'uint16');
          %remainder=65535;
          %for index=1:1:length(tempRawDataForCRC)
              %lookup=bitand((uint32(tempRawDataForCRC(index)^uint32(bitshift(remainder,(16-8),'uint32')))),255,'uint32');
              %remainder=uint32(crcTable(lookup)^bitshift(remainder,-8,'uint32'))
          %end
          %CRCcalc=uint32(bitand(remainder^0,65535,'uint32'));
          %if(CRC==CRCcalc)
              rawPackets=[rawPackets;tempPacket];
          %else
              %raise a flag
          %end
      end
  end
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%convert to floats here 
Uint16s=[[];[]];
temperature=[[];[]];
temperatureADC=[[];[]];
gimbalData=[[];[]];
successSensors=[];
for i=1:1:length(rawPackets)
    dataPointer=3;%exclude the length variable
    tempUint16=[];% stores the UINT16 variables 
    tempTemp=[];%temp array of temperature readings) 
    tempADCTemp=[];
    for j=1:1:N%N sensors //runs N times for sensors
        for k=1:1:9 
            tempUint16=[tempUint16 typecast(rawPackets(i,dataPointer:dataPointer+1),'UINT16')];
            dataPointer=dataPointer+2;
        end
        tempTemp=[tempTemp uint8(rawPackets(i,dataPointer))];
        dataPointer=dataPointer+1; 
        tempADCTemp=[tempADCTemp uint8(rawPackets(i,dataPointer))];%typecast(rawPackets(i,dataPointer:dataPointer+1),'UINT16')];
        dataPointer=dataPointer+1; 
          
    end
    tempGimbal16=[];
    if(gimbal==1)
        for j=1:1:4
            tempGimbal16=[tempGimbal16 typecast(rawPackets(i,dataPointer:dataPointer+3),'SINGLE')];%this should be floats 
            dataPointer=dataPointer+4;
        end
    end
    tempSuccess=0;
    if(successByte==1)
        tempSuccess=rawPackets(i,dataPointer);
        dataPointer=dataPointer+1;
        successSensors=[successSensors;tempSuccess];
    end
    gimbalData=[gimbalData;tempGimbal16];
    Uint16s=[Uint16s; tempUint16];
    temperature=[temperature; tempTemp];
    temperatureADC=[temperatureADC; tempADCTemp];
end

floats=[[];[]];
%twos comp
for i=1:1:length(Uint16s(:,1))
    tfloats=[];
    for j=1:1:length(Uint16s(1,:))%was a 2
        tfloats=[tfloats twoscomp(Uint16s(i,j))];%,Uint16s(i,j+1))];
    end
    floats=[floats;tfloats];
end

tempVal=[[];[]];%take out 
for i=1:1:length(temperature(:,1))
    ttemps=[];
    for j=1:1:length(temperature(1,:))  
        ttemps=[ttemps twoscompTemp(temperature(i))];
    end
    tempVal=[tempVal;ttemps];
end
% TODO scaling factor 
gyroScale=single(0.07);%this needs to change 
accScale=single(0.000786);
tempScale=1;
tempADCScale=1;

magScale1=single(0.0043);
magScale2=single(0.0049);
if(N>=1)
    gyroS1=floats(:,1:3)*gyroScale;
    magS1=floats(:,4:6);
    for i=1:1:length(magS1(:,1))
        magS1(i,:)=[magS1(i,1)*magScale1 magS1(i,2)*magScale1 magS1(i,3)*magScale2];
    end
    accS1=floats(:,7:9)*accScale;
    tempS1=tempVal(:,1)*tempScale;
    tempADCS1=temperatureADC(:,1)*tempADCScale;
end




if(gimbal==1)
    %get the gimbal data 
    rollEncoder=gimbalData(:,1);
    pitchEncoder=gimbalData(:,2);
    yawEncoder=gimbalData(:,3);
    yawVelocity=gimbalData(:,4);
    
end



end

