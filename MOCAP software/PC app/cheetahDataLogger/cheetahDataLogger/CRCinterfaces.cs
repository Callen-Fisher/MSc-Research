#region (c)Tellumat (Pty) Ltd 2009.
//
// Project              : TRACKS
// Classification       : RESTRICTED
// Filename             : CRCinterfaces.cs
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
    /// Interface classes for user defined CRC's
    /// </summary>
    public partial class CRC
    {
        /// <summary>
        /// mtek crc implementation for the interface
        /// </summary>
        public class CRC_mtek : CRCtype
        {
            public UInt32 polynomial { get { return 0x1021; } }
            public UInt32 initialRemainder { get { return 0; } }
            public UInt32 finalXORvalue { get { return 0; } }
            public bool reflectData { get { return true; } }
            public bool reflectRemainder { get { return true; } }
            public int dataSize { get { return 16; } }
        }

        /// <summary>
        /// cbads crc implementation for the interface
        /// </summary>
        public class CRC_cbads : CRCtype
        {
            public UInt32 polynomial { get { return 0x8005; } }
            public UInt32 initialRemainder { get { return 0; } }
            public UInt32 finalXORvalue { get { return 0; } }
            public bool reflectData { get { return false; } }
            public bool reflectRemainder { get { return false; } }
            public int dataSize { get { return 16; } }
        }

        /// <summary>
        /// data link crc implementation for the interface (CCITT)
        /// </summary>
        public class CRC_dataLink : CRCtype
        {
            public UInt32 polynomial { get { return 0x1021; } }
            public UInt32 initialRemainder { get { return 0xFFFF; } }
            public UInt32 finalXORvalue { get { return 0; } }
            public bool reflectData { get { return false; } }
            public bool reflectRemainder { get { return false; } }
            public int dataSize { get { return 16; } }
        }

        /// <summary>
        /// davis weather crc implementation for the interface
        /// </summary>
        public class CRC_Davis : CRCtype
        {
            public UInt32 polynomial { get { return 0x1021; } }
            public UInt32 initialRemainder { get { return 0; } }
            public UInt32 finalXORvalue { get { return 0; } }
            public bool reflectData { get { return false; } }
            public bool reflectRemainder { get { return false; } }
            public int dataSize { get { return 16; } }
        }

        /// <summary>
        /// Novatel crc implementation for the interface (CRC32)
        /// </summary>
        public class CRC_Novatel : CRCtype
        {
            public UInt32 polynomial { get { return 0x04C11DB7; } }
            public UInt32 initialRemainder { get { return 0; } }
            public UInt32 finalXORvalue { get { return 0; } }
            public bool reflectData { get { return true; } }
            public bool reflectRemainder { get { return true; } }
            public int dataSize { get { return 32; } }
        }

    }
}