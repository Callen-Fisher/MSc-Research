/*
 * SPIcomms.h
 *
 *  Created on: Nov 28, 2014
 *      Author: Callen Fisher
 */

#ifndef SPICOMMS_H_
#define SPICOMMS_H_
#include "stm32f4xx.h"
#include "stm32f4_discovery.h"
#include "hardwareSetup.h"
#include "FreeRTOS.h"
#include "vSensorTask.h"

void initSPI(void);
void setUpRegisters(void);
void chipSelect(void);
void chipDeSelect(void);
void setUpRegisters(void);
void delay100ms(void);
void delay130micro(void);
void chipEnable(void);
void chipDisable(void);
uint8_t writeSPIcomms(uint8_t regAdr, uint8_t data);
void readPacket(uint8_t* data,uint8_t length);
void setADDR(uint8_t Pipe, uint8_t one, uint8_t two, uint8_t three, uint8_t four, uint8_t five);
void setADDRshort(uint8_t Pipe, uint8_t one);
void sendRequest(uint8_t command);
void delay10micro(void);
void flushRX(void);
void flushTX(void);
void delay10ms(void);

uint8_t request;
uint8_t c;

#define PIPE0 0x0A
#define PIPE1 0x0B
#define PIPE2 0x0C
#define PIPE3 0x0D

#define PIPE1addr	0xB7
#define PIPE2addr	0xB8
#define PIPE3addr	0xB9

#define TXpipe 0x10

#endif /* SPICOMMS_H_ */
