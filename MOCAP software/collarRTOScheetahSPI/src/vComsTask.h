/*
 * vComsTask.h
 *
 *  Created on: Sep 29, 2014
 *      Author: Callen Fisher
 */

#ifndef VCOMSTASK_H_
#define VCOMSTASK_H_
#include "vSensorTask.h"
#include "task.h"

void vCommsTask( void *pvparameters );

xQueueHandle xCommsTransmitQueue;
xSemaphoreHandle xCommsTransmitSemaphore;

#endif /* VCOMSTASK_H_ */
