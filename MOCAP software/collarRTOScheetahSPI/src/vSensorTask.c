/*
 * vSensorTask.c
 *
 *  Created on: Sep 29, 2014
 *      Author: Callen Fisher
 */
#include "vSensorTask.h"


uint8_t adcVal;
int u=0;

void vSensorTask( void *pvparameters )
{
	vSensorTaskFrameReceived = xSemaphoreCreateMutex();
	xSemaphoreTake(vSensorTaskFrameReceived,0);

	//USART_Cmd(USART1,ENABLE);
	portTickType xLastWakeTime;
	//USART_Cmd(USART1,DISABLE);
	const portTickType xPeriod=LOOPTIME/portTICK_RATE_MS;
	xLastWakeTime=xTaskGetTickCount();
	for(;;)
	{
		vTaskDelayUntil(&xLastWakeTime,xPeriod);
		xSemaphoreTake(vSensorTaskFrameReceived, 1);


		//HACK DATA LOGGING LIGHT
		//STM_EVAL_LEDToggle(LED3);
		if(USART_GetFlagStatus(USART6,USART_FLAG_RXNE))//TOGGLE THE LOGGING LIGHT
		{
			GPIO_ToggleBits(GPIOE,GPIO_Pin_5);
			USART_ClearFlag(USART6,USART_FLAG_RXNE);
		}
		else
		{

		}
		//

		//success byte
		if(flagSuccessS1==1 && flagSuccessS2==1 && flagSuccessS3==1)
		{
			successByte=0b111;
		}
		else if(flagSuccessS1==0 && flagSuccessS2==1 && flagSuccessS3==1)
		{
			successByte=0b011;
		}
		else if(flagSuccessS1==1 && flagSuccessS2==0 && flagSuccessS3==1)
		{
			successByte=0b101;
		}
		else if(flagSuccessS1==1 && flagSuccessS2==1 && flagSuccessS3==0)
		{
			successByte=0b110;
		}
		else if(flagSuccessS1==0 && flagSuccessS2==0 && flagSuccessS3==1)
		{
			successByte=0b001;
		}
		else if(flagSuccessS1==1 && flagSuccessS2==0 && flagSuccessS3==0)
		{
			successByte=0b100;
		}
		else if(flagSuccessS1==0 && flagSuccessS2==1 && flagSuccessS3==0)
		{
			successByte=0b010;
		}
		else if(flagSuccessS1==0 && flagSuccessS2==0 && flagSuccessS3==0)
		{
			successByte=0b000;
		}

		flagSuccessS1=0;
		flagSuccessS2=0;
		flagSuccessS3=0;

		serialTerminal_packetize(dataS1,dataS2,dataS3,dataS4,sizeof(dataS1),sizeof(dataS2),sizeof(dataS3),sizeof(dataS4),successByte);
		for(q=0;q<sizeof(dataS1);q++)
		{
			dataS1[q]=0;
			dataS2[q]=0;
			dataS3[q]=0;
			dataS4[q]=0;
		}



		//TODO comms to the sensors:










		//read the data request
		s1_dma_init();//you dont need this
		USART_Cmd(USART1,ENABLE);//dont need
		USART_SendData(USART1,0x7E);//collar INEMO//dont need
		//S1
		setADDR(TXpipe,0xB3,0xB4,0xB5,0xB6,PIPE1addr);//set the address of the transmit pipe
		setADDR(PIPE0,0xB3,0xB4,0xB5,0xB6,PIPE1addr);//match that address for the ack packet, gets stored pipe 0
		sendRequest(readDataCommand);//sends data via SPI, gets stored in the buffer, readCommand gets replaced with ADC value
		chipEnable();//when enabled for 10micro starts sending over the air
		delay10micro();
		chipDisable();









		xSemaphoreTake(vSensorTaskFrameReceived, portMAX_DELAY);

		if(flagDataReceived==1)
		{
			readPacket(dataTemperature,sizeof(dataTemperature));
			STM_EVAL_LEDToggle(LED3);
			dataS1[18]=dataTemperature[0];
			dataS1[19]=dataTemperature[1];
		}
		flagDataReceived=0;
		flagRetriesFailed=0;
		flushTX();
		flushRX();

		//S2
		setADDR(TXpipe,0xB3,0xB4,0xB5,0xB6,PIPE2addr);
		setADDR(PIPE0,0xB3,0xB4,0xB5,0xB6,PIPE2addr);
		sendRequest(readDataCommand);
		chipEnable();
		delay10micro();
		chipDisable();

		xSemaphoreTake(vSensorTaskFrameReceived, portMAX_DELAY);

		if(flagDataReceived==1)
		{
			STM_EVAL_LEDToggle(LED3);
			readPacket(dataTemperature,sizeof(dataTemperature));
			dataS2[18]=dataTemperature[0];
			dataS2[19]=dataTemperature[1];
		}
		flagDataReceived=0;
		flagRetriesFailed=0;
		flushTX();
		flushRX();
		//S3
		setADDR(TXpipe,0xB3,0xB4,0xB5,0xB6,PIPE3addr);
		setADDR(PIPE0,0xB3,0xB4,0xB5,0xB6,PIPE3addr);
		sendRequest(readDataCommand);
		chipEnable();
		delay10micro();
		chipDisable();

		xSemaphoreTake(vSensorTaskFrameReceived, portMAX_DELAY);

		if(flagDataReceived==1)
		{
			STM_EVAL_LEDToggle(LED3);
			readPacket(dataTemperature,sizeof(dataTemperature));
			dataS3[18]=dataTemperature[0];
			dataS3[19]=dataTemperature[1];
		}
		flagDataReceived=0;
		flagRetriesFailed=0;
		flushTX();
		flushRX();

		//wait till data has been read
		xSemaphoreTake(vSensorTaskFrameReceived, 4/portTICK_RATE_MS);

		//request collar data
		//get the collar data
				DMA_Cmd(DMA2_Stream2,DISABLE);
				while (DMA_GetCmdStatus(DMA2_Stream2) != DISABLE)
				{
				}
				USART_Cmd(USART1,DISABLE);
				USART_DMACmd(USART1,USART_DMAReq_Rx,DISABLE);
				uint8_t s1_Rx_chars = RX_BUFFER_SIZE - DMA_GetCurrDataCounter(DMA2_Stream2);
				for(q=0;q<sizeof(Rx_Buffer);q++)
				{
					dataS4[q]=Rx_Buffer[q];
					Rx_Buffer[q]=0;
				}
				if(s1_Rx_chars==20)
					STM_EVAL_LEDOff(LED6);
				else
					STM_EVAL_LEDOn(LED6);
				//DMA_SetCurrDataCounter(DMA2_Stream2,19);

		//USART_Cmd(USART1,ENABLE);
		//USART_SendData(USART1,0xC2);

		//request S1 data

		setADDR(TXpipe,0xB3,0xB4,0xB5,0xB6,PIPE1addr);
		setADDR(PIPE0,0xB3,0xB4,0xB5,0xB6,PIPE1addr);
		sendRequest(requestSensorData);
		chipEnable();
		delay10micro();
		chipDisable();

		xSemaphoreTake(vSensorTaskFrameReceived, portMAX_DELAY);

		if(flagDataReceived==1)
		{
			//STM_EVAL_LEDToggle(LED6);
			readPacket(dataS1,sizeof(dataS1)-2);
			flagSuccessS1=1;
		}
		flagDataReceived=0;
		flagRetriesFailed=0;
		flushTX();
		flushRX();
		//request S2 data
		setADDR(TXpipe,0xB3,0xB4,0xB5,0xB6,PIPE2addr);
		setADDR(PIPE0,0xB3,0xB4,0xB5,0xB6,PIPE2addr);
		sendRequest(requestSensorData);
		chipEnable();
		delay10micro();
		chipDisable();

		xSemaphoreTake(vSensorTaskFrameReceived, portMAX_DELAY);

		if(flagDataReceived==1)
		{
			//STM_EVAL_LEDToggle(LED6);
			readPacket(dataS2,sizeof(dataS2)-2);
			flagSuccessS2=1;
		}
		flagDataReceived=0;
		flagRetriesFailed=0;
		flushTX();
		flushRX();
		//request S3 data
		setADDR(TXpipe,0xB3,0xB4,0xB5,0xB6,PIPE3addr);
		setADDR(PIPE0,0xB3,0xB4,0xB5,0xB6,PIPE3addr);
		sendRequest(requestSensorData);
		chipEnable();
		delay10micro();
		chipDisable();

		xSemaphoreTake(vSensorTaskFrameReceived, portMAX_DELAY);

		if(flagDataReceived==1)
		{
			//STM_EVAL_LEDToggle(LED6);
			readPacket(dataS3,sizeof(dataS3)-2);
			flagSuccessS3=1;
		}
		flagDataReceived=0;
		flagRetriesFailed=0;
		flushTX();
		flushRX();






		//CURRENT OF THE MOTOR HACK
		lengthS4=20;
		adcVal=adc_convert();//hack

		dataS4[19]=dataS4[18];
		dataS4[18]=adcVal;


	}
}
