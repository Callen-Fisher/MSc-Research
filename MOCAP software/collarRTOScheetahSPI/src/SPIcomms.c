/*
 * SPIcomms.c
 *
 *  Created on: Nov 28, 2014
 *      Author: Callen Fisher
 */

#include "SPIcomms.h"

void initSPI()//sets up the SPI between micro and the NRF chip
{
	GPIO_InitTypeDef GPIO_InitStructure;

	RCC_APB1PeriphClockCmd(RCC_APB1Periph_SPI2,ENABLE);
	RCC_AHB1PeriphClockCmd(RCC_AHB1Periph_GPIOD|RCC_AHB1Periph_GPIOB, ENABLE);
	RCC_APB2PeriphClockCmd(RCC_APB2Periph_SYSCFG,ENABLE);

	GPIO_PinAFConfig(GPIOB, GPIO_PinSource13, GPIO_AF_SPI2);//sck
	GPIO_PinAFConfig(GPIOB, GPIO_PinSource14, GPIO_AF_SPI2);//miso
	GPIO_PinAFConfig(GPIOB, GPIO_PinSource15, GPIO_AF_SPI2);//mosi

	GPIO_InitStructure.GPIO_Mode = GPIO_Mode_AF;
	GPIO_InitStructure.GPIO_Speed = GPIO_Speed_50MHz;
	GPIO_InitStructure.GPIO_OType = GPIO_OType_PP;
	GPIO_InitStructure.GPIO_PuPd  = GPIO_PuPd_DOWN;
	GPIO_InitStructure.GPIO_Pin = GPIO_Pin_13|GPIO_Pin_14|GPIO_Pin_15;
	GPIO_Init(GPIOB, &GPIO_InitStructure);

	SPI_InitTypeDef  SPI_InitStructure;

	SPI_I2S_DeInit(SPI2);
	SPI_InitStructure.SPI_Mode=SPI_Mode_Master;
	SPI_InitStructure.SPI_Direction = SPI_Direction_2Lines_FullDuplex;
	SPI_InitStructure.SPI_DataSize = SPI_DataSize_8b;
	SPI_InitStructure.SPI_CPOL = SPI_CPOL_Low;
	SPI_InitStructure.SPI_CPHA = SPI_CPHA_1Edge;
	SPI_InitStructure.SPI_NSS = SPI_NSS_Soft;
	SPI_InitStructure.SPI_BaudRatePrescaler = SPI_BaudRatePrescaler_2;
	SPI_InitStructure.SPI_FirstBit = SPI_FirstBit_MSB;
	SPI_InitStructure.SPI_CRCPolynomial = 7;

	SPI_Init(SPI2, &SPI_InitStructure);
	SPI_Cmd(SPI2, ENABLE);

	GPIO_InitTypeDef GPIO_InitStruct;
	GPIO_InitStruct.GPIO_Pin =GPIO_Pin_2|GPIO_Pin_1;//CSN =D2, CE=D1;
	GPIO_InitStruct.GPIO_Mode =GPIO_Mode_OUT;
	GPIO_InitStruct.GPIO_Speed =GPIO_Speed_50MHz;

	GPIO_Init(GPIOD, &GPIO_InitStruct);

	GPIO_InitStruct.GPIO_Pin =GPIO_Pin_0;//IRQ IRQ=B0
	GPIO_InitStruct.GPIO_Mode =GPIO_Mode_IN;
	GPIO_InitStruct.GPIO_Speed =GPIO_Speed_50MHz;
	GPIO_InitStruct.GPIO_PuPd=GPIO_PuPd_UP;
	GPIO_Init(GPIOD, &GPIO_InitStruct);

	SYSCFG_EXTILineConfig(EXTI_PortSourceGPIOD,EXTI_PinSource0);

	EXTI_InitTypeDef   EXTI_InitStructure;
	EXTI_InitStructure.EXTI_Line = EXTI_Line0;
	EXTI_InitStructure.EXTI_Mode = EXTI_Mode_Interrupt;
	EXTI_InitStructure.EXTI_Trigger = EXTI_Trigger_Falling;
	EXTI_InitStructure.EXTI_LineCmd = ENABLE;
	EXTI_Init(&EXTI_InitStructure);

	NVIC_InitTypeDef NVIC_InitStructure;
	NVIC_InitStructure.NVIC_IRQChannel = EXTI0_IRQn;
	NVIC_InitStructure.NVIC_IRQChannelPreemptionPriority = 5;
	NVIC_InitStructure.NVIC_IRQChannelSubPriority = 0x00;
	NVIC_InitStructure.NVIC_IRQChannelCmd = ENABLE;
	NVIC_Init(&NVIC_InitStructure);

	setUpRegisters();
}
void setUpRegisters()//sets up the nordic NRF chip to communicate as desired (parameters)
{
	chipDisable();
	delay100ms();
	//register 0x00
	//7 =0 only
	//6 =0=RX data ready interrupt
	//5 =0=TX data sent interrupt
	//4 =0=retries  interrupt
	//3 =1=enable CRC
	//2 =1=2 byte CRC
	//1 =1=power up
	//0 =1=PRX 	0=PTX (primary TX or RX)


	//when writing to a register:
	chipSelect();//first select the chip
	writeSPIcomms(0x00|0b00100000,0b00001110);//then use the SPI send function, need the reg address and data.
	chipDeSelect();//then you de-select the chip
	delay130micro();//delay before next operation
	delay130micro();




	//register 0x01
	//7 =0 only
	//6 =0 only
	//5 =0=disable ack pipe5
	//4 =0=disable ack pipe4
	//3 =0=disable ack pipe3
	//2 =0=disable ack pipe2
	//1 =0=disable ack pipe1
	//0 =1=enable ack pipe0
	chipSelect();
	writeSPIcomms(0x01|0b00100000,0b00000001);
	chipDeSelect();
	delay130micro();
	//register 0x02
	//7 =0 only
	//6 =0 only
	//5 =0=disable pipe5
	//4 =0=disable pipe4
	//3 =0=disable pipe3
	//2 =0=disable pipe2
	//1 =0=disable pipe1
	//0 =1=enable pipe0
	chipSelect();
	writeSPIcomms(0x02|0b00100000,0b00000001);
	chipDeSelect();
	delay130micro();
	//register 0x03
	//7 =0 only
	//6 =0 only
	//5 =0 only
	//4 =0 only
	//3 =0 only
	//2 =0 only
	//1 =1
	//0 =1 5 byte address
	chipSelect();
	writeSPIcomms(0x03|0b00100000,0b00000011);
	chipDeSelect();
	delay130micro();
	//register 0x04
	//7 =0
	//6 =1
	//5 =0
	//4 =1//auto retransmit wait 1500 us for all ACK payload sizes
	//3 =0
	//2 =0
	//1 =1
	//0 =1=3 retransmissions
	chipSelect();
	writeSPIcomms(0x04|0b00100000,0b01010011);
	chipDeSelect();
	delay130micro();
	//register 0x05
	//7 =0
	//6 =0
	//5 =0
	//4 =0
	//3 =0
	//2 =0
	//1 =1
	//0 =0 choose the frequency channel
	chipSelect();
	writeSPIcomms(0x04|0b00100000,0b00000010);
	chipDeSelect();
	delay130micro();
	//register 0x06
	//7 =0=allows continuous carrier transmit when high
	//6 =0 only
	//5 =0 only
	//4 =0 only
	//3 =1=2MBps
	//2 =1
	//1 =1=TX power of 0dBm
	//0 =1=dont care
	chipSelect();
	writeSPIcomms(0x06|0b00100000,0b00001111);
	chipDeSelect();
	delay130micro();


	//0x0A=p0 address(5 bytes)
	//0x0B=p1(5 bytes)
	//0x0C=p2(1 byte)
	//0x0D=p3(1 byte)
	//0x0E=p4(1 byte)
	//0x0F=p5(1 byte)

	//register 0x11
	//7 =0
	//6 =0
	//5 =0
	//4 =0
	//3 =0
	//2 =0
	//1 =0
	//0 =1
	chipSelect();
	writeSPIcomms(0x11|0b00100000,0b00000001);//1 byte P0
	chipDeSelect();
	delay130micro();
	//register 0x12
	//7 =0
	//6 =0
	//5 =0
	//4 =0
	//3 =0
	//2 =0
	//1 =0
	//0 =1
	//chipSelect();
	//writeSPIcomms(0x12|0b00100000,0b00000001);//1 byte P1
	//chipDeSelect();
	//delay130micro();

	//register 0x1C
	//7 =0 only
	//6 =0 only
	//5 =0
	//4 =0
	//3 =1 enable p3 dynamic payload length
	//2 =1 enable p2 dynamic payload length
	//1 =1 enable p1 dynamic payload length
	//0 =1 enable p0 dynamic payload length
	chipSelect();
	writeSPIcomms(0x1C|0b00100000,0b00000001);
	chipDeSelect();
	delay130micro();
	//register 0x1D
	//7 =0 only
	//6 =0 only
	//5 =0 only
	//4 =0 only
	//3 =0 only
	//2 =0 enable dynamic payload length
	//1 =1 enable payload with ACK
	//0 =0 disable the TX with no ACK command
	chipSelect();
	writeSPIcomms(0x1D|0b00100000,0b00000110);   ///////
	chipDeSelect();
	delay130micro();


	setADDR(TXpipe,0xB3,0xB4,0xB5,0xB6,PIPE1addr);
	setADDR(PIPE0,0xB3,0xB4,0xB5,0xB6,PIPE1addr);
}
void chipEnable()
{
	GPIO_SetBits(GPIOD, GPIO_Pin_1);
}
void chipDisable()
{
	GPIO_ResetBits(GPIOD, GPIO_Pin_1);
}
void chipSelect()
{
	GPIO_ResetBits(GPIOD, GPIO_Pin_2);
}
void chipDeSelect()
{
	GPIO_SetBits(GPIOD, GPIO_Pin_2);
}
void delay100ms()
{
	uint8_t p=0;
	uint8_t q=0;
	uint8_t r=0;
	uint8_t s=0;
	for(s=0;s<2;s++)
	{
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
	}
}
void delay10ms()
{
	uint8_t p=0;
	uint8_t q=0;
	uint8_t r=0;
	uint8_t s=0;
	for(s=0;s<1;s++)
	{
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
		for( p=0 ; p<254 ; p++)
		{
			for(r=0;r<254;r++)
			{
				q=p;
			}
		}
	}
}
void delay130micro()
{
	uint8_t p=0;
	uint8_t q=0;
	uint8_t r=0;
	uint8_t s=0;
	for(s=0;s<13;s++)
	{
	for( p=0 ; p<10 ; p++)
	{
		for(r=0;r<10;r++)
		{
			q=p;
		}
	}
	}
}
void delay10micro()
{
	uint8_t p=0;
	uint8_t q=0;
	uint8_t r=0;
	for( p=0 ; p<10 ; p++)
	{
		for(r=0;r<11;r++)
		{
			q=p;
		}
	}
}
uint8_t writeSPIcomms(uint8_t regAdr, uint8_t data)
{
	uint8_t dummyVar;
	int32_t val;

	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, regAdr);    //Sensor Address that we are WRITING to
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	dummyVar = SPI_I2S_ReceiveData(SPI2);

	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, data);    //Sensor Config
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	val = (u8)SPI_I2S_ReceiveData(SPI2);

	return (u8)val;
}
void setADDR(uint8_t Pipe, uint8_t one, uint8_t two, uint8_t three, uint8_t four, uint8_t five)
{
	//delay130micro();
	//delay10micro();
	chipSelect();
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, Pipe|0b00100000);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, one);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, two);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, three);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, four);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, five);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	chipDeSelect();
	//delay130micro();
	delay10micro();
}
void setADDRshort(uint8_t Pipe, uint8_t one)
{
	delay130micro();
	chipSelect();
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, Pipe|0b00100000);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, one);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	chipDeSelect();
	delay130micro();
}
void sendRequest(uint8_t command)
{
	chipSelect();
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, 0xA0);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, command);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	chipDeSelect();
	//delay130micro();
	delay10micro();
}
void readPacket(uint8_t* data,uint8_t length)
{
	//delay130micro();
	//delay10micro();
	chipSelect();
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, 0x61);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	uint8_t q;
	for(q=0;q<length;q++)
	{
		SPI_I2S_SendData(SPI2, 0x00);    //Sensor Address that we are WRITING to
		while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
		{
			/* Wait for data */
		}

	    data[q]=SPI_I2S_ReceiveData(SPI2);
		while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
		{
			/* Wait for all transmissions to complete */
		}
	}
	chipDeSelect();
	//delay130micro();
	delay10micro();
}
void flushTX()
{
	//delay130micro();
	//delay10micro();
	chipSelect();
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, 0xE1);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	chipDeSelect();
	//delay130micro();
	delay10micro();
}
void flushRX()
{
	//delay130micro();
	//delay10micro();
	chipSelect();
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	SPI_I2S_SendData(SPI2, 0xE2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_RXNE) == RESET)
	{
		/* Wait for data */
	}
	SPI_I2S_ReceiveData(SPI2);
	while (SPI_I2S_GetFlagStatus(SPI2, SPI_I2S_FLAG_TXE) == RESET)
	{
		/* Wait for all transmissions to complete */
	}
	chipDeSelect();
	//delay130micro();
	delay10micro();
}
