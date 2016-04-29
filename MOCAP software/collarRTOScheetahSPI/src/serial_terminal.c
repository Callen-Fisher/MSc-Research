/*
 * serial_terminal.c
 *
 *  Created on: Sep 15, 2014
 *      Author: Callen Fisher
 */

#include "serial_terminal.h"


static int8_t dePacketize(uint16_t entire_pkt_length, uint8_t* data_to_process, uint16_t* type_payload, uint8_t* extracted_data, uint16_t* data_length);

void serialTerminal_Init(void)
{
  initCRC();
}

void serialTerminal_packetize(uint8_t* payload_to_packS1,uint8_t* payload_to_packS2,uint8_t* payload_to_packS3,uint8_t* payload_to_packS4, uint16_t length_of_payloadS1,uint16_t length_of_payloadS2,uint16_t length_of_payloadS3,uint16_t length_of_payloadS4,uint8_t success)
{
  CommsTask_TransmitPacketStruct pkt_to_tx;
  uint16_t raw_data_pointer = 0;
  uint16_t packet_data_pointer = 4;
  uint16_t crcCalculated;
  uint8_t char_to_pack;
  uint8_t a=4;
  //put data
  for (raw_data_pointer = 0; raw_data_pointer < length_of_payloadS1; raw_data_pointer++, packet_data_pointer++,a++)
  { //for both bytes of the pkt_code,
    char_to_pack = payload_to_packS1[raw_data_pointer];
    t[a]=char_to_pack;
    if (char_to_pack == 0x7E)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5E;
    }
    else if (char_to_pack == 0x7D)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5D;
    }
    else
    {
      pkt_to_tx.data[packet_data_pointer] = char_to_pack;
    }
  }
  for (raw_data_pointer = 0; raw_data_pointer < length_of_payloadS2; raw_data_pointer++, packet_data_pointer++,a++)
  { //for both bytes of the pkt_code,
    char_to_pack = payload_to_packS2[raw_data_pointer];
    t[a]=char_to_pack;
    if (char_to_pack == 0x7E)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5E;
    }
    else if (char_to_pack == 0x7D)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5D;
    }
    else
    {
      pkt_to_tx.data[packet_data_pointer] = char_to_pack;
    }
  }
  for (raw_data_pointer = 0; raw_data_pointer < length_of_payloadS3; raw_data_pointer++, packet_data_pointer++,a++)
  { //for both bytes of the pkt_code,
    char_to_pack = payload_to_packS3[raw_data_pointer];
    t[a]=char_to_pack;
    if (char_to_pack == 0x7E)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5E;
    }
    else if (char_to_pack == 0x7D)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5D;
    }
    else
    {
      pkt_to_tx.data[packet_data_pointer] = char_to_pack;
    }
  }
  for (raw_data_pointer = 0; raw_data_pointer < length_of_payloadS4; raw_data_pointer++, packet_data_pointer++,a++)
  { //for both bytes of the pkt_code,
    char_to_pack = payload_to_packS4[raw_data_pointer];
    t[a]=char_to_pack;
    if (char_to_pack == 0x7E)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5E;
    }
    else if (char_to_pack == 0x7D)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5D;
    }
    else
    {
      pkt_to_tx.data[packet_data_pointer] = char_to_pack;
    }
  }



  t[a]=success;
  pkt_to_tx.data[packet_data_pointer]=success;//append which packets have been dropped and which haven't

  a++;
  packet_data_pointer++;



  uint16_t dataLen=a-4;//length of data



  union
  {
  	  char temp_char[2];
  	  uint16_t temp_16;
  }buffer_to_char_union1;
  buffer_to_char_union1.temp_16=0x0101;

  pkt_to_tx.data[0]=buffer_to_char_union1.temp_char[0];
  pkt_to_tx.data[1]=buffer_to_char_union1.temp_char[1];
  t[0]=buffer_to_char_union1.temp_char[0];
  t[1]=buffer_to_char_union1.temp_char[1];

  union
  {
     char temp_char[2];
     uint16_t temp_16;
  }buffer_to_char_union2;
  buffer_to_char_union2.temp_16=dataLen;

  pkt_to_tx.data[2]=buffer_to_char_union2.temp_char[0];
  pkt_to_tx.data[3]=buffer_to_char_union2.temp_char[1];
  t[2]=buffer_to_char_union2.temp_char[0];
  t[3]=buffer_to_char_union2.temp_char[1];
  //calculate CRC
  crcCalculated = (uint16_t) crcCalc(t, 0, a);//(uint16_t) crcCalc(pkt_to_tx.data, 4, packet_data_pointer-4);////////////////////COULD BE ERROR HERE

  //put CRC
  for (raw_data_pointer = 0; raw_data_pointer < 2; raw_data_pointer++, packet_data_pointer++)
  {
    char_to_pack = (uint8_t) ((crcCalculated >> 8 * (1 - raw_data_pointer)) & 0x00FF);
    if (char_to_pack == 0x7E)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5E;
    }
    else if (char_to_pack == 0x7D)
    {
      pkt_to_tx.data[packet_data_pointer++] = 0x7D;
      pkt_to_tx.data[packet_data_pointer] = 0x5D;
    }
    else
    {
      pkt_to_tx.data[packet_data_pointer] = char_to_pack;
    }
  }

  //set pkt length
  pkt_to_tx.bytes_to_tx = packet_data_pointer+2;
  //convert this to the TX buffer array
  CommsTask_TxBuffer_data[0]=0x7E;//start character

  //add the data (has the CRC,length and opcode in it)
  int i;
  for(i=1;i<pkt_to_tx.bytes_to_tx-1;i++)
  {
	  CommsTask_TxBuffer_data[i]=pkt_to_tx.data[i-1];
  }
  //add the terminating character
  CommsTask_TxBuffer_data[pkt_to_tx.bytes_to_tx-1]=0x7E;
  Tx_chars=pkt_to_tx.bytes_to_tx;


  //GPIO_ResetBits(GPIOE,GPIO_Pin_15);
  xSemaphoreGive(xCommsTransmitSemaphore);
  xQueueSend(xCommsTransmitQueue,CommsTask_TxBuffer_data,0);

  int w=0;
  for(w=0;w<sizeof(CommsTask_TxBuffer_data);w++)
  {
	  CommsTask_TxBuffer_data[w]=0;
  }
  //GPIO_SetBits(GPIOE,GPIO_Pin_15);
}

//Description: This function extracts the type, data and CRC payloads from the received packet,
//and performs a CRC check.
//Return: 0 for a sane packet, -1 for a broken one.
static int8_t dePacketize(uint16_t entire_pkt_length, uint8_t* data_to_process, uint16_t* type_payload, uint8_t* extracted_data, uint16_t* data_length)
{

//pkt structure:
//                 0        1  2     3   4        5,6,7...     n-3 n-2     n-1
//min length = start(1) + type(2) + length(2) + <data?>=length + crc(2) + end(1) = 8

  uint16_t crc_calculated;
  uint16_t crc_according_to_pkt;
  uint8_t raw_packet_pointer = 1; //start of payload (ie: the 1st byte of type)
  uint8_t extracted_payload_pointer = 0;

  if (entire_pkt_length < 8)//ie no data
  {
	  //STM_EVAL_LEDToggle(LED3);
    return -1;

  };

  if (data_to_process[0] != 0x7E)
  {
	  //STM_EVAL_LEDToggle(LED3);
    return -2;

  }

//7E -> 7D 5E
//7D -> 7D 5D

//extract payload:
  for (; raw_packet_pointer < entire_pkt_length - 1; extracted_payload_pointer++, raw_packet_pointer++)
  { //for all raw data characters
    if (data_to_process[raw_packet_pointer] == 0x7D)
    { //if the char is an escape char
      ++raw_packet_pointer; //look at the next received char
      if (data_to_process[raw_packet_pointer] == 0x5E)
      { //and if it is a 5E
        extracted_data[extracted_payload_pointer] = 0x7E; //replace it with a 7E
      }
      else if (data_to_process[raw_packet_pointer] == 0x5D)
      { //or if it is a 5D
        extracted_data[extracted_payload_pointer] = 0x7D; //replace it with a 7D
      }
      else
      { //or if it's something else,
        //printf("ESCAPE CHAR ERROR "); //we have an error
        //printf(" :: line %d in file %s \r\n", __LINE__ , __FILE__);

        return -3;

      }
    }
    else
    { //if not an escape char
      extracted_data[extracted_payload_pointer] = data_to_process[raw_packet_pointer]; //add normally.
      if(data_to_process[raw_packet_pointer]==0x7E)
      {
    	  //extracted_payload_pointer--;

    	  break;
      }
    }
  }


  *type_payload = extracted_data[0] | (extracted_data[1] << 8);

  *data_length = extracted_data[2] | (extracted_data[3] << 8);

  if (*data_length != extracted_payload_pointer - 6)//was 6 --changed to 5 due to -- in the above condition
  {

    if (*type_payload == 0x0101)
    {
     // *data_length -= 4;
    }
    else
    {
    	STM_EVAL_LEDToggle(LED5);
      return -4;


    }
  }

  crc_calculated = (uint16_t) crcCalc(extracted_data, 0, *data_length+4 );

  crc_according_to_pkt = extracted_data[extracted_payload_pointer - 1] | (extracted_data[extracted_payload_pointer - 2] << 8);//WAS -1 AND -2

  if (crc_calculated != crc_according_to_pkt)
  {
	  STM_EVAL_LEDToggle(LED6);
    return -5;

  }

  return 0;
}

void serial_terminal_ProcessPacketS1(uint16_t entire_pkt_length, uint8_t* data_to_process)
{
  uint16_t data_length = 0; //the length that the packet says it is
  uint16_t type_payload = 0;
  int8_t de_pack_return = 0;

  de_pack_return = dePacketize(entire_pkt_length, data_to_process, &type_payload, dataS1, &data_length);
  if (de_pack_return == 0)
  {

	int i;
    switch (type_payload)
    {
      case (0x0101):
    	//when depacetize left with the length and CRC and all that stuff, must remove it

    	for(i=4;i<data_length+4;i++)
    	{
    		dataS1[i-4]=dataS1[i];
    	}
        lengthS1=data_length;
        flagSuccessS1=1;
        break;
      default:
        break;
    }

  }
  else
  {
	  //STM_EVAL_LEDToggle(LED6);
	  flagSuccessS1=0;
	  flagSuccessS2=0;
	  flagSuccessS3=0;
  }
}
void serial_terminal_ProcessPacketS2(uint16_t entire_pkt_length, uint8_t* data_to_process)
{
  uint16_t data_length = 0; //the length that the packet says it is
  uint16_t type_payload = 0;
  int8_t de_pack_return = 0;

  de_pack_return = dePacketize(entire_pkt_length, data_to_process, &type_payload, dataS2, &data_length);

  if (de_pack_return == 0)
  {
	  int i;
	 // STM_EVAL_LEDToggle(LED3);

    switch (type_payload)
    {
      case (0x0101):

		for(i=4;i<data_length+4;i++)
		{
		    dataS2[i-4]=dataS2[i];
		}
        lengthS2=data_length;
        flagSuccessS2=1;
        break;
      default:
        break;
    }
  }
  else
  {
	  STM_EVAL_LEDToggle(LED6);
	  flagSuccessS1=0;
	  flagSuccessS2=0;
	  flagSuccessS3=0;

  }
}
void serial_terminal_ProcessPacketS3(uint16_t entire_pkt_length, uint8_t* data_to_process)
{
  uint16_t data_length = 0; //the length that the packet says it is
  uint16_t type_payload = 0;
  int8_t de_pack_return = 0;

  de_pack_return = dePacketize(entire_pkt_length, data_to_process, &type_payload, dataS3, &data_length);

  if (de_pack_return == 0)
  {
	  int i;

    switch (type_payload)
    {
      case (0x0101):

		for(i=4;i<data_length+4;i++)
		{
		    dataS3[i-4]=dataS3[i];
		}
        lengthS3=data_length;
        flagSuccessS3=1;

        break;
      default:
        break;
    }

  }
  else
  {
	  STM_EVAL_LEDToggle(LED6);
	  flagSuccessS1=0;
	  flagSuccessS2=0;
	  flagSuccessS3=0;
  }
}
