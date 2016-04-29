#region (c)Tellumat (Pty) Ltd 2009.
//
// Project              : TRACKS
// Classification       : RESTRICTED
// Filename             : CRC.cs
// Author               : C. van Aswegen

// The Copyright, manufacturing and patent rights stemming from this
// document in any form are vested in Tellumat (Pty) Ltd. Tellumat has
// ceded these rights to its clients where contractually agreed.
// 
#endregion

using System;
using System.Collections.Generic;
using System.Text;

namespace WindowsFormsApplication1
{
    /// <summary>
    /// A universal CRC library
    /// </summary>
    public partial class CRC
    {
        #region INTERFACES
        /// <summary>
        /// Interface class for used by the constructor to set up the crc class
        /// </summary>
        public interface CRCtype
        {
            UInt32 polynomial { get; }
            UInt32 initialRemainder { get; }
            UInt32 finalXORvalue { get; }
            bool reflectData { get; }
            bool reflectRemainder { get; }
            int dataSize { get; }
        }
        #endregion

        #region FIELDS
        UInt32[] crcTable;
        int mask;
        #endregion

        #region PROPERTIES
        /// <summary>
        /// Current ongoing CRC value as the bytes are added. Initial value is 0x0000
        /// </summary>
        UInt32 _CRCvalue;
        public UInt32 CRCvalue 
        {
            get
            {
                if (tCRC.reflectRemainder)
                {
                    return ((UInt32)reflect(_CRCvalue, tCRC.dataSize) ^ tCRC.finalXORvalue);
                }
                else
                {
                    return (_CRCvalue ^ tCRC.finalXORvalue);
                }
            }
            set
            {
                _CRCvalue = value;
            }
        } 

        public CRCtype tCRC { get; private set; }
        #endregion

        #region CONSTRUCTOR
        /// <summary>
        /// Constructor sets up the look up table based on the interface class implemented
        /// </summary>
        /// <param name="ctype"></param>
        public CRC(CRCtype ctype)
        {
            tCRC = ctype;
            _CRCvalue = ctype.initialRemainder;
            crcTable = new UInt32[256];
            UInt32 remainder;
            UInt32 dividend;
            int bit;
            int topbit;

            mask = math.intpow(2, ctype.dataSize) - 1;
            topbit = 1 << (ctype.dataSize - 1);
            // Compute the remainder of each possible dividend.
            for (dividend = 0; dividend < 256; ++dividend)
            {
                remainder = (UInt32)((dividend << (ctype.dataSize - 8)) & mask);
                for (bit = 8; bit > 0; --bit)
                {
                    // Try to divide the current data bit.
                    if ((remainder & topbit) != 0)
                    {
                        remainder = (UInt32)(((remainder << 1) ^ ctype.polynomial) & mask);
                    }
                    else
                    {
                        remainder = (UInt32)((remainder << 1) & mask);
                    }
                }
                crcTable[dividend] = remainder;
            }
        }
        #endregion

        #region PUBLIC METHODS
        /// <summary>
        /// Calculate CRC of input buffer.
        /// </summary>
        /// <param name="buffer">byte array</param>
        /// <param name="Length">length of byte array</param>
        /// <returns>the calculated CRC</returns>
        public UInt32 CalcCrc(byte[] buffer, int Length)
        {
            int remainder = (int)tCRC.initialRemainder;
            int lookup;
            int i;
            // Divide the message by the polynomial, a byte at a time.
            for (i = 0; i < Length; ++i)//USE TO START AT 0
            {
                if (tCRC.reflectData)
                {
                    lookup = (int)(reflect((UInt32)buffer[i], 8) ^ (remainder >> (tCRC.dataSize - 8))) & 0xFF;
                }
                else
                {
                    lookup = (buffer[i] ^ (remainder >> (tCRC.dataSize - 8))) & 0xFF;
                }
                remainder = (int)((crcTable[lookup] ^ (remainder << 8)) & mask);
            }
            // The final remainder is the CRC.
            if (tCRC.reflectRemainder)
            {
                return (reflect((UInt32)remainder, tCRC.dataSize) ^ tCRC.finalXORvalue);
            }
            else
            {
                UInt32 val = (UInt32)(remainder ^ tCRC.finalXORvalue);
                return val;
            }
        }

        /// <summary>
        /// Calculate CRC of input string
        /// </summary>
        /// <param name="buffer">string</param>
        /// <returns>the calculated CRC</returns>
        public UInt32 CalcCrc(string buffer)
        {

            return CalcCrc(ASCIIEncoding.UTF8.GetBytes(buffer), buffer.Length);
        }
        ///// <summary>
        ///// Adds the byte to CRCvalue
        ///// </summary>
        ///// <param name="b">byte to be added to CRCvalue</param>
        ///// <returns>Return the current value of CCRCvalue after adding incoming byte</returns>
        //public int AddByte(byte b)
        //{
        //    int data;
        //    if (tCRC.reflectData)
        //    {
        //        data = (reflect(b, 8) ^ (_CRCvalue >> (tCRC.dataSize - 8))) & 0xFF;
        //    }
        //    else
        //    {
        //        data = (b ^ (_CRCvalue >> (tCRC.dataSize - 8))) & 0xFF;
        //    }
        //    _CRCvalue = (crcTable[data] ^ (_CRCvalue << 8)) & mask;
        //    if (tCRC.reflectRemainder)
        //    {
        //        return (reflect(_CRCvalue, tCRC.dataSize) ^ tCRC.finalXORvalue);
        //    }
        //    else
        //    {
        //        return (_CRCvalue ^ tCRC.finalXORvalue);
        //    }
        //}

        /// <summary>
        /// Resets the CRC
        /// </summary>
        public void Clear()
        {
            _CRCvalue = tCRC.initialRemainder;
        }

        #endregion 

        #region PRIVATE METHODS
        /// <summary>
        /// Does a bitwise reflection on data
        /// </summary>
        /// <param name="data">data to reflect</param>
        /// <param name="nBits">data type size</param>
        /// <returns>reflected value</returns>
        UInt32 reflect(UInt32 data, int nBits)
        {
            int reflection = 0;
            int bit;
            // Reflect the data about the center bit.
            for (bit = 0; bit < nBits; ++bit)
            {
                // If the LSB bit is set, set the reflection of it.
                if ((data & 0x01) > 0)
                {
                    reflection |= (1 << ((nBits - 1) - bit));
                }
                data = (data >> 1);
            }
            return (UInt32)reflection;
        }
        #endregion
    }
}
