/*
 * vComsTask.c
 *
 *  Created on: Sep 29, 2014
 *      Author: Callen Fisher
 */

#include "vComsTask.h"
int p=0;
int count=0;

void vCommsTask(void *pvparameters)
{
	//uint16_t length;
	xCommsTransmitQueue = xQueueCreate(12, sizeof(CommsTask_TxBuffer));
	xCommsTransmitSemaphore = xSemaphoreCreateMutex();

	comms_dma_init();

	for(;;)
	{
		xSemaphoreTake(xCommsTransmitSemaphore,portMAX_DELAY);
		STM_EVAL_LEDToggle(LED4);
		xQueueReceive(xCommsTransmitQueue, CommsTask_TxBuffer,0);

		count=0;
		for(p=0;p<sizeof(CommsTask_TxBuffer);p++)
		{
			if(CommsTask_TxBuffer[p]==0x7E)count++;
			if(count>2)CommsTask_TxBuffer[p]=0;
		}
		DMA_Cmd(PCDMASTREAM,DISABLE);
		while(DMA_GetCmdStatus(PCDMASTREAM) != DISABLE)
		{

		}

		DMA_SetCurrDataCounter(PCDMASTREAM,sizeof(CommsTask_TxBuffer));//length);
		DMA_Cmd(PCDMASTREAM,ENABLE);
		while(DMA_GetCmdStatus(PCDMASTREAM) != ENABLE)
		{

		}
		xSemaphoreTake(xCommsTransmitSemaphore, portMAX_DELAY);
		//clear the buffer
		for(p=0;p<sizeof(CommsTask_TxBuffer);p++)
		{
			CommsTask_TxBuffer[p]=0;
		}

	}
}
