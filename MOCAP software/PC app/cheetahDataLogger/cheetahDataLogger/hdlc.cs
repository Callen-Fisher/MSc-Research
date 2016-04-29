#region (c)Tellumat (Pty) Ltd 2009.
//
// Project              : TRACKS
// Classification       : RESTRICTED
// Filename             : HDLC.cs
// Author               : C. van Aswegen

// The Copyright, manufacturing and patent rights stemming from this
// document in any form are vested in Tellumat (Pty) Ltd. Tellumat has
// ceded these rights to its clients where contractually agreed.
// 
#endregion

#region PREPROCESSOR
#define DEBUG_HDLC
#endregion

using System;
using System.Collections.Generic;
using System.Text;
using System.Diagnostics;


namespace WindowsFormsApplication1
{
    /// <summary>
    /// Processes HDLC frame data and builds valid data packets to be processed
    /// </summary>
    public class HDLC
    {
        #region FIELDS
        int crcFailCount = 0;
        int escapeCharFail = 0;
        bool RxSOF;
        bool ESCFlag;
        List<byte> RxBuffer;
        List<byte> originalFrameReceived;
        int ByteCount;
        #endregion

        #region CONSTRUCTOR
        /// <summary>
        /// Creates the HDLC class
        /// </summary>
        public HDLC()
        {
            RxBuffer = new List<byte>();
            originalFrameReceived = new List<byte>();
            RxSOF = false;
            ESCFlag = false;
            ByteCount = 0;
        }
        #endregion

        #region PUBLIC METHODS
        /// <summary>
        /// Build HDLC frame data packet of the data
        /// </summary>
        /// <param name="data">data packet to frame</param>
        /// <returns>framed data</returns>
        public byte[] BuildTxFrame(byte[] data)
        {
            if (data == null)
            {
                return null;
            }
            List<byte> TXdata = new List<byte>();
            CRC.CRC_dataLink Txlink = new CRC.CRC_dataLink();
            CRC TxCRC = new CRC(Txlink);
            uint TxCRCval = TxCRC.CalcCrc(data, data.Length);
            //SOF
            TXdata.Add(0x7E);
            //Add data to TXbuffer
            for (int i = 0; i < data.Length; i++)
            {
                AddToTxBuffer(ref TXdata, data[i]);
            }
            //Add CRC to TXbuffer
            byte HighCRC = (byte)((TxCRCval >> 8) & 0xFF);
            byte LowCRC = (byte)(TxCRCval & 0xFF);
            AddToTxBuffer(ref TXdata, HighCRC);
            AddToTxBuffer(ref TXdata, LowCRC);
            //EOF
            TXdata.Add(0x7E);
            return TXdata.ToArray();
        }

        /// <summary>
        /// Build up the RX data on a byte by byte basis. Return byte array if complete frame is received, otherwise null.
        /// </summary>
        /// <param name="InByte">input data</param>
        /// <returns>valid data, else null</returns>
        public byte[] ProcessRxByte(byte InByte, byte[] originalFrame)
        {
            originalFrameReceived.Add(InByte);
            ByteCount++;
            
            if (InByte == 0x7E)
            {
                originalFrame[0] = InByte;
                originalFrame[ByteCount] = InByte;
                if (RxSOF)
                {
                    if (ByteCount == 1)
                    {
                        //caught the packet mid frame
                        ByteCount = 0;
                        #if (DEBUG && DEBUG_HDLC)
                            Debug.WriteLine(DateTime.Now.ToString() + " - HDLC - Double SOF detected");
                        #endif
                    }
                    else
                    {
                        CRC RxCRC = new CRC(new CRC.CRC_dataLink());

                       // ##############################################################################
                        if (RxCRC.CalcCrc(RxBuffer.ToArray(), RxBuffer.Count)==0) 
                        {
                            RxSOF = false;
                            //RxBuffer.RemoveRange(RxBuffer.Count - 2, 2);//???????????????????????????????????????
                            return RxBuffer.ToArray();
                        }
                        else
                        {
                            RxSOF = false;
                            #if (DEBUG && DEBUG_HDLC)
                                Debug.WriteLine(DateTime.Now.ToString() + " - HDLC - Invalid CRC");
                            #endif
                            RxBuffer.Clear();
                            crcFailCount++;
                            //throw new InvalidCRCException("CRC value incorrect");
                           //return RxBuffer.ToArray();/////THIS MUST COME OUT
                        }
                    }
                }
                else
                {
                    RxSOF = true;
                    RxBuffer.Clear();
                    ESCFlag = false;
                    ByteCount = 0;
                   // return RxBuffer.ToArray();/////THIS MUST COME OUT
                    return null;
                }
            }
            else
            {
                originalFrame[ByteCount] = InByte;
                if (RxSOF == false)
                {
                    return null;
                }
                if (InByte == 0x7D)
                {
                    ESCFlag = true;
                    return null;
                }
                if (ESCFlag)
                {
                    if (InByte == 0x5E)
                    {
                        RxBuffer.Add(0x7E);
                        ESCFlag = false;
                        return null;
                    }
                    else if (InByte == 0x5D)
                    {
                        RxBuffer.Add(0x7D);
                        ESCFlag = false;
                        return null;
                    }
                    else
                    {
                        RxSOF = false;
                        #if (DEBUG && DEBUG_HDLC)
                            Debug.WriteLine(DateTime.Now.ToString() + " - HDLC - Incorrect ESC char - 0x" + InByte.ToString("X2"));
                        #endif
                        escapeCharFail++;
                        //throw new InvalidESCException("Incorrect ESC sequence received");
                    }
                }
                else
                {
                    RxBuffer.Add(InByte);
                }
            }
            return null;
        }
        #endregion

        public int getNumberOfCRCFails()
        {
            return crcFailCount;
        }

        public int getNumberOfescapeCharFails()
        {
            return escapeCharFail;
        }

        public void resetCounters()
        {
            crcFailCount = 0;
            escapeCharFail = 0;
        }

        #region PRIVATE METHODS
        /// <summary>
        /// Add to the TX buffer and does escaping if neccessary
        /// </summary>
        /// <param name="TXdata">ref TXdata buffer</param>
        /// <param name="data">byte to add</param>
        private void AddToTxBuffer(ref List<byte> TXdata, byte data)
        {
            if (data == 0x7E)
            {
                //ESC character
                TXdata.Add(0x7D);
                TXdata.Add(0x5E);
            }
            else if (data == 0x7D)
            {
                //ESC character
                TXdata.Add(0x7D);
                TXdata.Add(0x5D);
            }
            else
            {
                TXdata.Add(data);
            }
        }
        #endregion
    }

    #region EXCEPTION CLASSES
    /// <summary>
    /// Exception thrown if the packet fails CRC
    /// </summary>
    public class InvalidCRCException : Exception
    {
        public InvalidCRCException() { }
        public InvalidCRCException(string message) : base(message) { }
        public InvalidCRCException(string message, System.Exception inner) : base(message, inner) { }
    }
    /// <summary>
    /// Exception thrown if the incorrect escape sequence is received
    /// </summary>
    public class InvalidESCException : Exception
    {
        public InvalidESCException() { }
        public InvalidESCException(string message) : base(message) { }
        public InvalidESCException(string message, System.Exception inner) : base(message, inner) { }
    }
    #endregion
}
