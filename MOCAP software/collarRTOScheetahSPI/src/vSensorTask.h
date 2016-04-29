/*
 * vSensorTask.h
 *
 *  Created on: Sep 29, 2014
 *      Author: Callen Fisher
 */

#ifndef VSENSORTASK_H_
#define VSENSORTASK_H_
#include "FreeRTOS.h"
#include "stm32f4xx.h"
#include "queue.h"
#include "semphr.h"
#include "serial_terminal.h"
#include "hardwareSetup.h"
#include "adcCurrent.h"
#include "SPIcomms.h"
#include "XBEE.h"

#define LengthInemoPacket     20
#define LOOPTIME    		  10
uint8_t  successByte;

int q;

#define readDataCommand		  0x7E
#define requestSensorData	  0xC3

uint8_t flagDataReceived;
uint8_t flagDataSent;
uint8_t flagRetriesFailed;

void vSensorTask( void *pvparameters );

xSemaphoreHandle vSensorTaskFrameReceived;

#endif /* VSENSORTASK_H_ */
