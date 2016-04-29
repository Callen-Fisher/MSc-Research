/*
 * hardwareSetup.h
 *
 *  Created on: Sep 15, 2014
 *      Author: Callen Fisher
 */

#ifndef HARDWARESETUP_H_
#define HARDWARESETUP_H_
#include "stm32f4xx.h"
#include "stm32f4_discovery.h"

// The # defines:
#define RX_BUFFER_SIZE 				20
#define TX_BUFFER_SIZE 				120
//COMMS TO PC
#define DMA2BUS						RCC_AHB1Periph_DMA2
#define	PCUSARTBUS					RCC_APB2Periph_USART6
#define PCPORT						GPIOC
#define PCPINSOURCE					GPIO_PinSource6
#define PCPINSOURCERX				GPIO_PinSource7
#define PCALTERNATEFUNCTION			GPIO_AF_USART6
#define PCTXPIN						GPIO_Pin_6
#define PCRXPIN						GPIO_Pin_7
#define	PCBAUDRATE					230400
#define PCUSART						USART6
#define PORTABUS					RCC_AHB1Periph_GPIOA
#define PORTBBUS					RCC_AHB1Periph_GPIOB
#define PORTCBUS					RCC_AHB1Periph_GPIOC
#define PORTDBUS  					RCC_AHB1Periph_GPIOD

//PCDMA
#define PCDMASTREAM					DMA2_Stream6
#define PCDMACHANNEL				DMA_Channel_5
#define PCDMABASEADDRESS			USART6_BASE
#define PCDMABUFFER					(uint32_t) CommsTask_TxBuffer
#define PCDMABUFFERSIZE				sizeof(CommsTask_TxBuffer)



//DATA ARRAYS
uint8_t CommsTask_TxBuffer[TX_BUFFER_SIZE];
uint8_t CommsTask_TxBuffer_data[TX_BUFFER_SIZE];
uint8_t t[TX_BUFFER_SIZE];

uint8_t dataS1[RX_BUFFER_SIZE];
uint8_t dataS2[RX_BUFFER_SIZE];
uint8_t dataS3[RX_BUFFER_SIZE];
uint8_t dataS4[RX_BUFFER_SIZE];
uint8_t dataTemperature[2];
//CHARACTERS RECEIVED VARIABLES
uint8_t Tx_chars;


//FUNCTIONS
void comms_dma_init(void);
void setupCommsHardware();
#endif
