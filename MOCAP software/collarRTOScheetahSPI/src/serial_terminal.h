/*
 * serial_terminal.h
 *
 *  Created on: Sep 15, 2014
 *      Author: Callen Fisher
 */

#ifndef SERIAL_TERMINAL_H_
#define SERIAL_TERMINAL_H_

#include "CRC.h"
#include "hardwareSetup.h"
#include "stm32f4xx.h"
#include "string.h"
#include "vComsTask.h"

typedef struct {
  uint8_t data[TX_BUFFER_SIZE - 2];

  uint16_t bytes_to_tx;
} CommsTask_TransmitPacketStruct;

//FLAGS FOR SUCCESSFUL PACKETS
uint8_t flagSuccessS1;
uint8_t flagSuccessS2;
uint8_t flagSuccessS3;

//DATA RECEIVED LENGTHS
uint16_t lengthS1;
uint16_t lengthS2;
uint16_t lengthS3;
uint16_t lengthS4;

void serialTerminal_Init();
void serialTerminal_packetize(uint8_t* payload_to_packS1,uint8_t* payload_to_packS2,uint8_t* payload_to_packS3,uint8_t* payload_to_packS4, uint16_t length_of_payloadS1, uint16_t length_of_payloadS2, uint16_t length_of_payloadS3, uint16_t length_of_payloadS4,uint8_t success);
void SerialTerminal_BeginTx(void);

#endif /* SERIAL_TERMINAL_H_ */
