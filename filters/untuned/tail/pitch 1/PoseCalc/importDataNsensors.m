

N=4;%the number of sensors that were data logged (gimbal logs 1 and app logs 4)
gimbal=0;
successByte=1;% is there a success byte appended to the packet (tells which sensors were dropped

if(N==4)%if N=4 then using the collar to log so there will be a success byte appended 
    successByte=1;
end

if(gimbal==1)% gimbal always data loggs one sensor
    N=1;
    successByte=0;
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
accScale=single(0.000786/6);
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
if(N>=2)
    gyroS2=floats(:,10:12)*gyroScale;
    magS2=floats(:,13:15);
    for i=1:1:length(magS2(:,1))
        magS2(i,:)=[magS2(i,1)*magScale1 magS2(i,2)*magScale1 magS2(i,3)*magScale2];
    end
    accS2=floats(:,16:18)*accScale;
    tempS2=tempVal(:,2)*tempScale;
    tempADCS2=temperatureADC(:,2)*tempADCScale;
end
if(N>=3)
    gyroS3=floats(:,19:21)*gyroScale;
    magS3=floats(:,22:24);
    for i=1:1:length(magS3(:,1))
        magS3(i,:)=[magS3(i,1)*magScale1 magS3(i,2)*magScale1 magS3(i,3)*magScale2];
    end
    accS3=floats(:,25:27)*accScale;
    tempS3=tempVal(:,3)*tempScale;
    tempADCS3=temperatureADC(:,3)*tempADCScale;
end
if(N>=4)
    gyroS4=floats(:,28:30)*gyroScale;
    magS4=floats(:,31:33);
    for i=1:1:length(magS4(:,1))
        magS4(i,:)=[magS4(i,1)*magScale1 magS4(i,2)*magScale1 magS4(i,3)*magScale2];
    end
    accS4=floats(:,34:36)*accScale;
    tempS4=tempVal(:,4)*tempScale;
    tempADCS4=temperatureADC(:,4)*tempADCScale;
end




if(gimbal==1)
    %get the gimbal data 
    rollEncoder=gimbalData(:,1);
    pitchEncoder=gimbalData(:,2);
    yawEncoder=gimbalData(:,3);
    yawVelocity=gimbalData(:,4);
    
end


% working=0;//rather turn off the update part of the filter 
% for i=2:1:length(successSensors)-1
%     if(successSensors(i)==7)
%         %do nothing
%     elseif(successSensors(i)==6)
%         %sensor 3 is not working
%         
%         %find when next the sensor is working:
%         for j=i+1:1:length(successSensors)-1
%             if(successSensors(j)==6)
%                 working=j;
%                 break;
%             end
%         end
%         
%         
%         if(working~=0)
%             aveAccS3=(accS3(:,i-1)+accS3(:,working))/2;
%             aveMagS3=(magS3(:,i-1)+magS3(:,working))/2;
%             aveGyroS3=(gyroS3(:,i-1)+gyroS3(:,working))/2;
%             aveTempS3=(tempS3(i-1)+tempS3(working));
%             aveTempADCS3=(tempADCS3(i-1)+tempADCS3(working));
%             
%             for j=i:1:working-1
%                 accS3(:,j)=aveAccS3(:);
%                 magS3(:,j)=aveMagS3(:);
%                 gyroS3(:,j)=aveGyroS3(:);
%                 tempS3(j)=aveTempS3;
%                 tempADCS3(j)=aveTempADCS3(:);
%                 
%                 successByte(j)=7;
%             end
%         end 
%         working=0;
%     elseif(successSensors(i)==5)
%         %sensor 2 is not working
%     elseif(successSensors(i)==4)
%         %sensor 2 and 3 are not working
%     elseif(successSensor(i)==3)
%         %sensor 1 is not working
%     elseif(successSensor(i)==2)
%         %sensor 1 and 3 are not working
%     elseif(successSensor(i)==1)
%         %sensor 1 and 2 are not working 
%     end
% end





clear tempVal ttemps i j floats tfloats temperature dataPointer tempUint16 
clear gimbalData N Uint16s ans fid gimbal inByte k opcode rawPackets successByte 
clear successSensor temp tempGimbal16 tempPacket tempSuccess tempTemp tempRawDataForCRC
clear tempADCTemp temperatureADC working









