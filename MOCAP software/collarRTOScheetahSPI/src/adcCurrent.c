/*
 * adcCurrent.c
 *
 *  Created on: Nov 24, 2014
 *      Author: Callen Fisher
 */
#include "adcCurrent.h"

void setupADC()
{//hack for the current sensor

	//PB0
	//adc12_IN8
	ADC_InitTypeDef ADC_init_structure;
	GPIO_InitTypeDef GPIO_initStructre;
	RCC_APB2PeriphClockCmd(RCC_APB2Periph_ADC2,ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1ENR_GPIOBEN,ENABLE);

	GPIO_initStructre.GPIO_Pin = GPIO_Pin_0;
	GPIO_initStructre.GPIO_Mode = GPIO_Mode_AN;
	GPIO_initStructre.GPIO_PuPd = GPIO_PuPd_NOPULL;
	GPIO_Init(GPIOB,&GPIO_initStructre);

	ADC_DeInit();
	ADC_init_structure.ADC_DataAlign = ADC_DataAlign_Right;
	ADC_init_structure.ADC_Resolution = ADC_Resolution_8b;
	ADC_init_structure.ADC_ContinuousConvMode = ENABLE;
	ADC_init_structure.ADC_ExternalTrigConv = ADC_ExternalTrigConv_T1_CC1;
	ADC_init_structure.ADC_ExternalTrigConvEdge = ADC_ExternalTrigConvEdge_None;
	ADC_init_structure.ADC_NbrOfConversion = 1;
	ADC_init_structure.ADC_ScanConvMode = DISABLE;
	ADC_Init(ADC2,&ADC_init_structure);
	ADC_Cmd(ADC2,ENABLE);
	ADC_RegularChannelConfig(ADC2,ADC_Channel_8,1,ADC_SampleTime_144Cycles);
}
uint8_t adc_convert(){
 ADC_SoftwareStartConv(ADC2);
 while(!ADC_GetFlagStatus(ADC2, ADC_FLAG_EOC));
 return (uint8_t)(ADC_GetConversionValue(ADC2));
}
