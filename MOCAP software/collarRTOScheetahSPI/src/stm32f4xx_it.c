/**
*****************************************************************************
**
**  File        : stm32f4xx_it.c
**
**  Abstract    : Main Interrupt Service Routines.
**                This file provides template for all exceptions handler and
**                peripherals interrupt service routine.
**
**  Environment : Atollic TrueSTUDIO(R)
**                STMicroelectronics STM32F4xx Standard Peripherals Library
**
**  Distribution: The file is distributed “as is,” without any warranty
**                of any kind.
**
**  (c)Copyright Atollic AB.
**  You may use this file as-is or modify it according to the needs of your
**  project. Distribution of this file (unmodified or modified) is not
**  permitted. Atollic AB permit registered Atollic TrueSTUDIO(R) users the
**  rights to distribute the assembled, compiled & linked contents of this
**  file as part of an application binary file, provided that it is built
**  using the Atollic TrueSTUDIO(R) toolchain.
**
**
*****************************************************************************
*/

/* Includes ------------------------------------------------------------------*/
#include "stm32f4xx_it.h"
#include "FreeRTOS.h"
#include "task.h"
#include "vComsTask.h"
#include "vSensorTask.h"
#include "SPIcomms.h"
/* Private typedef -----------------------------------------------------------*/
/* Private define ------------------------------------------------------------*/
/* Private macro -------------------------------------------------------------*/
/* Private variables ---------------------------------------------------------*/
/* Private function prototypes -----------------------------------------------*/
/* Private functions ---------------------------------------------------------*/

/******************************************************************************/
/*            Cortex-M4 Processor Exceptions Handlers                         */
/******************************************************************************/
// ----------------------------------------------------------------------------
void EXTI0_IRQHandler(void)//IRQ line
{
	portBASE_TYPE woken = pdFALSE;
	if(EXTI_GetITStatus(EXTI_Line0))
	{
		chipSelect();
		uint8_t val=writeSPIcomms(0x07,0x07);//check the interrupt
		chipDeSelect();
		//delay130micro();
		if((val&0b01000000)==0b01000000)
		{
			//data received
			flagDataReceived=1;

			chipSelect();
			writeSPIcomms(0x07|0b00100000,0b01110000);
			chipDeSelect();
			//delay130micro();
			EXTI_ClearITPendingBit(EXTI_Line0);

			xSemaphoreGiveFromISR(vSensorTaskFrameReceived, &woken);
			portEND_SWITCHING_ISR(woken);
		}
		else if((val&0b00100000)==0b00100000)
		{
			//transmit complete-data sent

			flagDataSent=1;
			chipSelect();
			writeSPIcomms(0x07|0b00100000,0b01110000);
			chipDeSelect();
			//delay130micro();
			EXTI_ClearITPendingBit(EXTI_Line0);

			xSemaphoreGiveFromISR(vSensorTaskFrameReceived, &woken);
			portEND_SWITCHING_ISR(woken);
		}
		else if((val&0b00010000)==0b00010000)
		{
			flagRetriesFailed=1;
			//number of retries failed
			STM_EVAL_LEDToggle(LED5);
			chipSelect();
			writeSPIcomms(0x07|0b00100000,0b01110000);
			chipDeSelect();
			//delay130micro();
			EXTI_ClearITPendingBit(EXTI_Line0);

			xSemaphoreGiveFromISR(vSensorTaskFrameReceived, &woken);
			portEND_SWITCHING_ISR(woken);
		}

	}
}

void DMA2_Stream6_IRQHandler(void)
{
  portBASE_TYPE higher_task_woken;

  //Check to see that the IRQ was caused by DMA1_Stream_5 Transfer_Complete
  if (DMA_GetITStatus(PCDMASTREAM, DMA_IT_TCIF6) != RESET) {
    xSemaphoreGiveFromISR(xCommsTransmitSemaphore, &higher_task_woken);
    DMA_ClearITPendingBit(PCDMASTREAM, DMA_IT_TCIF6);
    //portEND_SWITCHING_ISR(higher_task_woken);
  }
}

/**
  * @brief   This function handles NMI exception.
  * @param  None
  * @retval None
  */
void NMI_Handler(void)
{
}

/**
  * @brief  This function handles Hard Fault exception.
  * @param  None
  * @retval None
  */
void HardFault_Handler(void)
{
  /* Go to infinite loop when Hard Fault exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles Memory Manage exception.
  * @param  None
  * @retval None
  */
void MemManage_Handler(void)
{
  /* Go to infinite loop when Memory Manage exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles Bus Fault exception.
  * @param  None
  * @retval None
  */
void BusFault_Handler(void)
{
  /* Go to infinite loop when Bus Fault exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles Usage Fault exception.
  * @param  None
  * @retval None
  */
void UsageFault_Handler(void)
{
  /* Go to infinite loop when Usage Fault exception occurs */
  while (1)
  {
  }
}

/**
  * @brief  This function handles SVCall exception.
  * @param  None
  * @retval None
  */
//void SVC_Handler(void)
//{
//}

/**
  * @brief  This function handles Debug Monitor exception.
  * @param  None
  * @retval None
  */
void DebugMon_Handler(void)
{
}


// JEK -> The following two function are handled by FreeRTOS.  See line 225 
// in port.c inside of FreeRTOS.

/**
  * @brief  This function handles PendSVC exception.
  * @param  None
  * @retval None
  */
//void PendSV_Handler(void)
//{
//}

/**
  * @brief  This function handles SysTick Handler.
  * @param  None
  * @retval None
  */
//void SysTick_Handler(void)
//{
//}

/******************************************************************************/
/*                 STM32F4xx Peripherals Interrupt Handlers                   */
/*  Add here the Interrupt Handler for the used peripheral(s) (PPP), for the  */
/*  available peripheral interrupt handler's name please refer to the startup */
/*  file (startup_stm32f4xx.s).                                               */
/******************************************************************************/

/**
  * @brief  This function handles PPP interrupt request.
  * @param  None
  * @retval None
  */
/*void PPP_IRQHandler(void)
{
}*/

