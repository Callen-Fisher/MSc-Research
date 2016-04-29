/*
 * XBEE.h
 *
 *  Created on: Mar 25, 2015
 *      Author: User
 */

#ifndef XBEE_H_
#define XBEE_H_
#include "stm32f4xx.h"

void setupCollarInemo(void);
void s1_dma_init(void);
uint8_t Rx_Buffer[19];
#endif /* XBEE_H_ */
