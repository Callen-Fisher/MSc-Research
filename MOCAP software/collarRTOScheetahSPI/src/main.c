
#include "FreeRTOS.h"
#include "task.h"
#include "queue.h"
#include "semphr.h"
#include "stm32f4xx.h"
#include "stm32f4_discovery.h"
#include "utils.h"
#include "hardwareSetup.h"
#include "vComsTask.h"
#include "vSensorTask.h"
#include "main.h"
#include "SPIcomms.h"
#include "serial_terminal.h"
#include "XBEE.h"

int main( void )
{

	delay100ms();

	STM_EVAL_LEDInit(LED3);
	STM_EVAL_LEDInit(LED4);
	STM_EVAL_LEDInit(LED5);
	STM_EVAL_LEDInit(LED6);
	STM_EVAL_LEDToggle(LED3);


	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOE,ENABLE);
	GPIO_InitTypeDef GPIO_InitStruct;
	GPIO_InitStruct.GPIO_Pin =GPIO_Pin_5;
	GPIO_InitStruct.GPIO_Mode =GPIO_Mode_OUT;
	GPIO_InitStruct.GPIO_Speed =GPIO_Speed_100MHz;
	GPIO_Init(GPIOE, &GPIO_InitStruct);

	initSPI();
	delay100ms();
	EXTI_ClearITPendingBit(EXTI_Line0);
	flushRX();
	flushTX();

	setupCommsHardware();
	serialTerminal_Init();
	setupADC();
	setupCollarInemo();
	xTaskCreate( vSensorTask, ( signed char * ) "SENSOR_TASK", 1200, NULL,3, NULL );//800//1000
	xTaskCreate( vCommsTask, ( signed char * ) "COMMS_TASK", 912, NULL,2, NULL );//512//712
    vTaskStartScheduler();

    for( ;; );  
}
void vApplicationTickHook( void )
{
    ++tickTime;
}
void vApplicationIdleHook( void )
{
    ++u64IdleTicksCnt;
}
void vApplicationMallocFailedHook( void )
{
    configASSERT( 0 );
}
