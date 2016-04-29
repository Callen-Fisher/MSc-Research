/*
 * XBEE.c
 *
 *  Created on: Mar 25, 2015
 *      Author: User
 */
#include "XBEE.h"

void setupCollarInemo()
{
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2, ENABLE);
	RCC_APB2PeriphClockCmd(RCC_APB2Periph_USART1,ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOB,ENABLE);
	GPIO_InitTypeDef GPIO_InitStruct;

	USART_InitTypeDef USART_InitStructure;

	GPIO_InitStruct.GPIO_Pin =GPIO_Pin_7;
	GPIO_InitStruct.GPIO_Mode =GPIO_Mode_AF;
	GPIO_InitStruct.GPIO_OType = GPIO_OType_PP;
	GPIO_InitStruct.GPIO_Speed =GPIO_Speed_100MHz;
	GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_UP;
	GPIO_Init(GPIOB, &GPIO_InitStruct);
	GPIO_PinAFConfig(GPIOB, GPIO_PinSource7, GPIO_AF_USART1);


	GPIO_InitStruct.GPIO_Pin = GPIO_Pin_6;
	GPIO_InitStruct.GPIO_Mode = GPIO_Mode_AF;
	GPIO_InitStruct.GPIO_OType = GPIO_OType_PP;
	GPIO_InitStruct.GPIO_Speed = GPIO_Speed_100MHz;
	GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_UP;
	GPIO_Init(GPIOB, &GPIO_InitStruct);
	GPIO_PinAFConfig(GPIOB, GPIO_PinSource6, GPIO_AF_USART1);

	USART_InitStructure.USART_BaudRate=230400;
	USART_InitStructure.USART_WordLength=USART_WordLength_8b;
	USART_InitStructure.USART_StopBits=USART_StopBits_1;
	USART_InitStructure.USART_Parity=USART_Parity_No;
	USART_InitStructure.USART_Mode=USART_Mode_Rx|USART_Mode_Tx;
	USART_InitStructure.USART_HardwareFlowControl=USART_HardwareFlowControl_None;

	USART_Init(USART1, &USART_InitStructure);

	USART_Cmd(USART1,DISABLE);
}
void s1_dma_init()//receive data from sensor
{
	//setup the DMA
	DMA_InitTypeDef DMA_InitStructure;

	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_DMA2, ENABLE);

	DMA_DeInit(DMA2_Stream2);
	while (DMA_GetCmdStatus(DMA2_Stream2) != DISABLE)
	{
	}

	//Initialize the structure to default values
	DMA_StructInit(&DMA_InitStructure);
	DMA_InitStructure.DMA_Channel = DMA_Channel_4;
	DMA_InitStructure.DMA_PeripheralBaseAddr = USART1_BASE + (uint32_t)0x04; //USART Data Register
	DMA_InitStructure.DMA_Memory0BaseAddr = (uint32_t)Rx_Buffer;
	DMA_InitStructure.DMA_DIR = DMA_DIR_PeripheralToMemory;
	DMA_InitStructure.DMA_BufferSize = sizeof(Rx_Buffer);
	DMA_InitStructure.DMA_PeripheralInc = DMA_PeripheralInc_Disable;
	DMA_InitStructure.DMA_MemoryInc = DMA_MemoryInc_Enable;
	DMA_InitStructure.DMA_PeripheralDataSize = DMA_PeripheralDataSize_Byte;
	DMA_InitStructure.DMA_MemoryDataSize = DMA_MemoryDataSize_Byte;
	DMA_InitStructure.DMA_Mode = DMA_Mode_Circular;
	DMA_InitStructure.DMA_Priority = DMA_Priority_VeryHigh;
	DMA_InitStructure.DMA_FIFOMode = DMA_FIFOMode_Disable;
	DMA_InitStructure.DMA_FIFOThreshold = DMA_FIFOThreshold_Full;
	DMA_InitStructure.DMA_MemoryBurst = DMA_MemoryBurst_Single;
	DMA_InitStructure.DMA_PeripheralBurst = DMA_PeripheralBurst_Single;
	DMA_Init(DMA2_Stream2, &DMA_InitStructure);





	DMA_Cmd(DMA2_Stream2, ENABLE);

		while (DMA_GetCmdStatus(DMA2_Stream2) != ENABLE)
		{
		}
		USART_DMACmd(USART1, USART_DMAReq_Rx, ENABLE);
}
void DMA2_Stream2_IRQHandler(void) // USART3 RX DMA Transfer Complete Interrupt
{
  if (DMA_GetITStatus(DMA2_Stream2, DMA_IT_TC))
  {
	// ----------

    DMA_ClearITPendingBit(DMA2_Stream2, DMA_IT_TC);
  }
}
