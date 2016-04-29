/*
 * hardwareSetup.c
 *
 *  Created on: Sep 15, 2014
 *      Author: Callen Fisher
 */
#include "hardwareSetup.h"


void setupCommsHardware()
{
	RCC_AHB1PeriphClockCmd(DMA2BUS, ENABLE);
	RCC_APB2PeriphClockCmd(PCUSARTBUS,ENABLE);
	RCC_AHB1PeriphClockCmd(PORTABUS|PORTBBUS|PORTCBUS|PORTDBUS,ENABLE);

	//COMMS TO THE PC
	GPIO_InitTypeDef GPIO_InitStruct;
	USART_InitTypeDef USART_InitStructure;

	GPIO_PinAFConfig(PCPORT, PCPINSOURCE, PCALTERNATEFUNCTION);

	GPIO_PinAFConfig(PCPORT, PCPINSOURCERX, PCALTERNATEFUNCTION);
	GPIO_InitStruct.GPIO_Pin =PCTXPIN|PCRXPIN;
	GPIO_InitStruct.GPIO_Mode =GPIO_Mode_AF;
	GPIO_InitStruct.GPIO_OType = GPIO_OType_PP;
	GPIO_InitStruct.GPIO_Speed =GPIO_Speed_100MHz;
	GPIO_InitStruct.GPIO_PuPd = GPIO_PuPd_UP;
	GPIO_Init(PCPORT, &GPIO_InitStruct);

	USART_StructInit(&USART_InitStructure);
	USART_InitStructure.USART_BaudRate=PCBAUDRATE;
	USART_InitStructure.USART_WordLength=USART_WordLength_8b;
	USART_InitStructure.USART_StopBits=USART_StopBits_1;
	USART_InitStructure.USART_Parity=USART_Parity_No;
	USART_InitStructure.USART_Mode=USART_Mode_Tx|USART_Mode_Rx;
	USART_InitStructure.USART_HardwareFlowControl=USART_HardwareFlowControl_None;

	USART_Init(PCUSART, &USART_InitStructure);
	USART_Cmd(PCUSART,ENABLE);

}
void comms_dma_init()//transmit data via XBEE
{
	DMA_InitTypeDef DMA_InitStructure;
	NVIC_InitTypeDef NVIC_InitStructure;
	RCC_AHB1PeriphClockCmd(DMA2BUS, ENABLE);

	DMA_DeInit(PCDMASTREAM);
	while (DMA_GetCmdStatus(PCDMASTREAM) != DISABLE)
	{
	}


    DMA_StructInit(&DMA_InitStructure);
    DMA_InitStructure.DMA_Channel = PCDMACHANNEL;
    DMA_InitStructure.DMA_PeripheralBaseAddr = PCDMABASEADDRESS + 0x04; //USART Data Register
    DMA_InitStructure.DMA_Memory0BaseAddr = PCDMABUFFER;
    DMA_InitStructure.DMA_DIR = DMA_DIR_MemoryToPeripheral;
    DMA_InitStructure.DMA_BufferSize = PCDMABUFFERSIZE;
    DMA_InitStructure.DMA_PeripheralInc = DMA_PeripheralInc_Disable;
    DMA_InitStructure.DMA_MemoryInc = DMA_MemoryInc_Enable;
    DMA_InitStructure.DMA_PeripheralDataSize = DMA_PeripheralDataSize_Byte;
    DMA_InitStructure.DMA_MemoryDataSize = DMA_MemoryDataSize_Byte;
    DMA_InitStructure.DMA_Mode = DMA_Mode_Normal;
    DMA_InitStructure.DMA_Priority = DMA_Priority_Low;
    DMA_InitStructure.DMA_FIFOMode = DMA_FIFOMode_Disable;
    DMA_InitStructure.DMA_FIFOThreshold = DMA_FIFOThreshold_Full;
    DMA_InitStructure.DMA_MemoryBurst = DMA_MemoryBurst_Single;
    DMA_InitStructure.DMA_PeripheralBurst = DMA_PeripheralBurst_Single;
    DMA_Init(PCDMASTREAM, &DMA_InitStructure);

    DMA_ITConfig(PCDMASTREAM, DMA_IT_TC, ENABLE);

    //Configure DMA TX Stream interrupt in the NVIC
    NVIC_InitStructure.NVIC_IRQChannel = DMA2_Stream6_IRQn;
    NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = 5;
    NVIC_InitStructure.NVIC_IRQChannelSubPriority = 0x00;
    NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
    NVIC_Init(&NVIC_InitStructure);


    USART_DMACmd(PCUSART, USART_DMAReq_Tx, ENABLE);

    DMA_Cmd(PCDMASTREAM,ENABLE);
    while (DMA_GetCmdStatus(PCDMASTREAM) != ENABLE)
    {
    }

}
