/*
 * adcCurrent.h
 *
 *  Created on: Nov 24, 2014
 *      Author: Callen Fisher
 */

#ifndef ADCCURRENT_H_
#define ADCCURRENT_H_
#include "FreeRTOS.h"
#include "stm32f4xx.h"
#include "stdint.h"

void setupADC(void);
uint8_t adc_convert(void);

#endif /* ADCCURRENT_H_ */
