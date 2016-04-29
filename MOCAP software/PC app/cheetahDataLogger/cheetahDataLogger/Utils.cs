#region (c)Tellumat (Pty) Ltd 2009.
//
// Project              : TRACKS
// Classification       : RESTRICTED
// Filename             : Utils.cs
// Author               : C. van Aswegen

// The Copyright, manufacturing and patent rights stemming from this
// document in any form are vested in Tellumat (Pty) Ltd. Tellumat has
// ceded these rights to its clients where contractually agreed.
// 
#endregion

using System;
using System.Collections.Generic;
using System.Text;
using System.IO;
using System.Diagnostics;
using System.Runtime.Serialization.Formatters.Binary;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Windows.Forms;
using System.Drawing;
using System.Data;

namespace WindowsFormsApplication1
{
    /// <summary>
    /// Static methods for the conversion of CTime
    /// </summary>
    public class Ctime
    {
        /// <summary>
        /// converts ctime (secs since 1/1/1970 12:00am) to DateTime
        /// </summary>
        /// <param name="CTime">ctime value</param>
        /// <returns>DateTime value</returns>
        public static DateTime CTimeToDateTime(uint CTime)
        {
            TimeSpan span = TimeSpan.FromTicks(CTime * TimeSpan.TicksPerSecond);
            DateTime t = new DateTime(1970, 1, 1).Add(span);
            return t;
        }

        /// <summary>
        /// converts DateTime to ctime (secs since 1/1/1970 12:00am)
        /// </summary>
        /// <param name="dt">DateTime value</param>
        /// <returns>ctime value</returns>
        public static uint DateTimeToCTime(DateTime dt)
        {
            DateTime t = new DateTime(1970, 1, 1);
            TimeSpan span = new TimeSpan(dt.Ticks - t.Ticks);
            return ((uint)(span.Ticks/TimeSpan.TicksPerSecond));
        }
    }

    //public class Thread
    //{
    //    public static void InvokeControlAction<t>(t cont, Action<t> action) where t : Control
    //    {
    //        if (cont.InvokeRequired)
    //        {
    //            cont.Invoke(new Action<t, Action<t>>(InvokeControlAction),
    //                      new object[] { cont, action });
    //        }
    //        else
    //        { action(cont); }
    //    }
    //}

    /// <summary>
    /// Static methods for math utilities
    /// </summary>
    public class math
    {
        /// <summary>
        /// Calculates the power of an integer a^b
        /// </summary>
        /// <param name="a">base</param>
        /// <param name="b">power</param>
        /// <returns>calculated value</returns>
        public static int intpow(int a, int b)
        {
            int result = 1;
            while (b != 0)
            {
                if ((b & 1) != 0)
                    result *= a;
                a *= a;
                b >>= 1;
            }
            return result;
        }

        /// <summary>
        /// Convert radians to degrees
        /// </summary>
        /// <param name="rad">radians</param>
        /// <returns>degrees</returns>
        public static double Rad2Deg(double rad)
        {
            return rad * (180 / Math.PI);
        }

        /// <summary>
        /// Convert degrees to radians
        /// </summary>
        /// <param name="deg">degrees</param>
        /// <returns>radians</returns>
        public static double Deg2Rad(double deg)
        {
            return deg * (Math.PI / 180);
        }

        /// <summary>
        /// Convert Milli-arcseconds to degrees
        /// </summary>
        /// <param name="Val">milli-arcseconds</param>
        /// <returns>degrees</returns>
        public static double MilliArcSecs2Deg(Int32 Val)
        {
            return (double)Val / 3600000;
        }
        /// <summary>
        /// Convert degrees to Milli-arcseconds
        /// </summary>
        /// <param name="Val">degrees</param>
        /// <returns>milli-arcseconds</returns>
        public static Int32 Deg2MilliArcSecs(double Val)
        {
            return (Int32)(Val * 3600000);
        }

        /// <summary>
        /// Convert radians to Milli-arcseconds
        /// </summary>
        /// <param name="Val">radians</param>
        /// <returns>milli-arcseconds</returns>
        public static Int32 Rad2MilliArcSecs(double rad)
        {
            return (Int32)((rad * (180 / Math.PI)) * 3600000);
        }

        /// <summary>
        /// Convert north and east velocity to resultant velocity
        /// </summary>
        /// <param name="VelocityNorth">north velocity component</param>
        /// <param name="VelocityEast">east velocity component</param>
        /// <returns>resultant velocity</returns>
        public static double Velocity2D(double VelocityNorth, double VelocityEast)
        {
            return Math.Sqrt((VelocityNorth * VelocityNorth) +
                    (VelocityEast * VelocityEast));
        }

        /// <summary>
        /// Convert north and east velocity to heading
        /// </summary>
        /// <param name="VelocityNorth">north velocity component</param>
        /// <param name="VelocityEast">east velocity component</param>
        /// <returns>resultant heading</returns>
        public static double Heading(double VelocityNorth, double VelocityEast)
        {
            if (VelocityEast == 0)
            {
                if (VelocityNorth >= 0)
                    return 0;
                else
                    return 180;
            }
            double VectAngle = Rad2Deg(Math.Atan(Math.Abs(VelocityNorth) / Math.Abs(VelocityEast)));
            if ((VelocityEast > 0) && (VelocityNorth >= 0))
            {
                return 90 - VectAngle;
            }
            else if ((VelocityEast > 0) && (VelocityNorth <= 0))
            {
                return 90 + VectAngle;
            }
            else if ((VelocityEast < 0) && (VelocityNorth <= 0))
            {
                return 270 - VectAngle;
            }
            else
            {
                return 270 + VectAngle;
            }
            //return (Rad2Deg(Math.Atan2(VelocityNorth, VelocityEast)) * -1) - 90;
        }

        /// <summary>
        /// Calculates the distance between 2 points
        /// </summary>
        /// <param name="lat1">latitude in radians</param>
        /// <param name="lon1">longitude in radians</param>
        /// <param name="lat2">latitude in radians</param>
        /// <param name="lon2">longitude in radians</param>
        /// <returns>distance in meters</returns>
        public static double DistanceBetweenWaypointsRad(double lat1, double lon1, double lat2, double lon2)
        {
            double SinVal = Math.Sin(lat1) * Math.Sin(lat2);
            double CosVal = Math.Cos(lat1) * Math.Cos(lat2) * Math.Cos((-1 * lon2) - (-1 * lon1));
            double CosDist = Math.Acos((CosVal) + (SinVal));
            return CosDist * 6378137;
        }

        /// <summary>
        /// Calculates the distance between 2 points
        /// </summary>
        /// <param name="lat1">latitude in degrees</param>
        /// <param name="lon1">longitude in degrees</param>
        /// <param name="lat2">latitude in degrees</param>
        /// <param name="lon2">longitude in degrees</param>
        /// <returns>distance in meters</returns>
        public static double DistanceBetweenWaypointsDeg(double lat1, double lon1, double lat2, double lon2)
        {
            return DistanceBetweenWaypointsRad(Deg2Rad(lat1), Deg2Rad(lon1), Deg2Rad(lat2), Deg2Rad(lon2));
        }

        /// <summary>
        /// Calculates the bearing between 2 points
        /// </summary>
        /// <param name="lat1">latitude in degrees</param>
        /// <param name="lon1">longitude in degrees</param>
        /// <param name="lat2">latitude in degrees</param>
        /// <param name="lon2">longitude in degrees</param>
        /// <returns>bearing in degrees</returns>
        public static double BearingBetweenWaypoints(double lat1, double lon1, double lat2, double lon2)
        {
            double bearing = AzimuthBetweenWaypoints(lat1, lon1, lat2, lon2);
            if (bearing < 0)
            {
                return 360 + bearing;
            }
            else
            {
                return bearing;
            }
        }

        /// <summary>
        /// Calculates the azimuth between 2 points
        /// </summary>
        /// <param name="lat1">latitude in degrees</param>
        /// <param name="lon1">longitude in degrees</param>
        /// <param name="lat2">latitude in degrees</param>
        /// <param name="lon2">longitude in degrees</param>
        /// <returns>azimuth in degrees (-180 to +180)</returns>
        public static double AzimuthBetweenWaypoints(double lat1, double lon1, double lat2, double lon2)
        {
            double Azimuth = Rad2Deg(Math.Atan2(Math.Sin(Deg2Rad(lon2) - Deg2Rad(lon1)) * Math.Cos(Deg2Rad(lat2)),
                                Math.Cos(Deg2Rad(lat1)) * Math.Sin(Deg2Rad(lat2)) - Math.Sin(Deg2Rad(lat1)) * Math.Cos(Deg2Rad(lat2)) * Math.Cos(Deg2Rad(lon2) - Deg2Rad(lon1))));
            return Azimuth;
        }

        /// <summary>
        /// Calculates the elevation between 2 points
        /// </summary>
        /// <param name="lat1">latitude in degrees</param>
        /// <param name="lon1">longitude in degrees</param>
        /// <param name="h1">altitude in meters</param>
        /// <param name="lat2">latitude in degrees</param>
        /// <param name="lon2">longitude in degrees</param>
        /// <param name="h2">altitude in meters</param>
        /// <returns>elevation in degrees</returns>
        public static double ElevationBetweenCoordinates(double lat1, double lon1, double h1, double lat2, double lon2, double h2)
        {
            double R = 6378137;
            double d = DistanceBetweenWaypointsDeg(lat1, lon1, lat2, lon2);
            double hg1 = h1 + R;
            double hg2 = h2 + R;
            double md = Math.Sqrt((hg1 * hg1) + (d * d));
            double dt = hg2 - md;
            return Rad2Deg(Math.Atan2(dt, d));
        }
        /// <summary>
        /// Calculates the distance between 2 point (pythagoras)
        /// </summary>
        /// <param name="X1">Point 1 - X</param>
        /// <param name="Y1">Point 1 - Y</param>
        /// <param name="X2">Point 2 - X</param>
        /// <param name="Y2">Point 2 - Y</param>
        /// <returns>distance</returns>
        public static double DistanceBetweenPoints(double X1, double Y1, double X2, double Y2)
        {
            return Math.Sqrt(((X2 - X1) * (X2 - X1)) + ((Y2 - Y1) * (Y2 - Y1)));
        }

        /// <summary>
        /// Calculates the distance between 2 points (pythagoras)
        /// </summary>
        /// <param name="pos1">first point</param>
        /// <param name="pos2">second point</param>
        /// <returns>distance</returns>
        public static double DistanceBetweenPoints(PointF pos1, PointF pos2)
        {
            return Math.Sqrt(((pos2.X - pos1.X) * (pos2.X - pos1.X)) + ((pos2.Y - pos1.Y) * (pos2.Y - pos1.Y)));
        }

        /// <summary>
        /// return the point on a circle based on the angle away from the X axis
        /// </summary>
        /// <param name="X">circle center</param>
        /// <param name="Y">circle center</param>
        /// <param name="radius">circle radius</param>
        /// <param name="angle">angle of the point</param>
        /// <returns>point of the circle for the angle position</returns>
        public static PointF PointOnCircle(double X, double Y, double radius, double angle)
        {
            float pX, pY;
            pX = (float)(X + (radius * Math.Cos(Deg2Rad(angle))));
            pY = (float)(Y + (radius * Math.Sin(Deg2Rad(angle))));
            return new PointF(pX, pY);
        }

        /// <summary>
        /// Converts meters per second to knots
        /// </summary>
        /// <param name="speed"></param>
        /// <returns></returns>
        public static double MetersPerSecond2Knots(double speed)
        {
            return speed * 1.943844;
        }

        /// <summary>
        /// converts knots to meters per second
        /// </summary>
        /// <param name="speed">speed in knots</param>
        /// <returns>speed in meters per second</returns>
        public static double Knots2MetersPerSecond(double speed)
        {
            return speed * 0.514444;
        }

        /// <summary>
        /// Converts meters per second to feet per second
        /// </summary>
        /// <param name="speed">speed in meters per second</param>
        /// <returns>speed in feet per second</returns>
        public static double MetersPerSecond2FeetPerSecond(double speed)
        {
            return speed * 3.280840;
        }

        /// <summary>
        /// Converts feet per second to meters per second
        /// </summary>
        /// <param name="speed">speed in feet per second</param>
        /// <returns>speed in meters per second</returns>
        public static double FeetPerSecond2MetersPerSecond(double speed)
        {
            return speed * 0.3048;
        }

        /// <summary>
        /// Convert meters per second per feet per minute
        /// </summary>
        /// <param name="speed">speed in meters per second</param>
        /// <returns>speed in feet per minute</returns>
        public static double MetersPerSecond2FeetPerMinute(double speed)
        {
            return speed * 3.280840 * 60;
        }

        /// <summary>
        /// Convert miles per hour to meters per second
        /// </summary>
        /// <param name="speed">speed in miles per hour</param>
        /// <returns>speed in metes per second</returns>
        public static double MilesPerHour2MetersPerSecond(double speed)
        {
            return speed * 0.44704;
        }
        /// <summary>
        /// Convert miles per hour to kilometers per hour
        /// </summary>
        /// <param name="speed">speed in miles per hour</param>
        /// <returns>speed in kilometers per hour</returns>
        public static double MilesPerHour2KilometersPerHour(double speed)
        {
            return speed * 1.609344;
        }
        /// <summary>
        /// Convert miles per hour to knots
        /// </summary>
        /// <param name="speed">speed in miles per hour</param>
        /// <returns>speed in knots</returns>
        public static double MilesPerHour2Knots(double speed)
        {
            return speed * 0.868976;
        }
        /// <summary>
        /// Convert degrees fahrenheit to degrees celcius
        /// </summary>
        /// <param name="fahrenheit">temperature in degrees fahrenheit</param>
        /// <returns>temperature in degrees celcius</returns>
        public static double Fahrenheit2Celcius(double fahrenheit)
        {
            return (fahrenheit - 32) * ((double)5 / (double)9);
        }
        /// <summary>
        /// Convert pressure in inches mercury to millibar
        /// </summary>
        /// <param name="pressure">pressure in inches mercury</param>
        /// <returns>pressure in millibar</returns>
        public static double inHg2mbar(double pressure)
        {
            return pressure * 33.85;
        }
    }

    /// <summary>
    /// Static methods for working with individual digits in a number
    /// </summary>
    public class Numbers
    {
        /// <summary>
        /// Gets the digis that make up the whole part of the number
        /// </summary>
        /// <param name="val">number to get the digits</param>
        /// <param name="MinDigits">minimum number of digits to return</param>
        /// <returns>array of digits form right to left from the decimal point</returns>
        public static int[] WholeDigits(float val, int MinDigits)
        {
            int[] wholes = WholeDigits(val);
            if (wholes.Length <= MinDigits)
            {
                int[] digitsVals = new int[MinDigits];
                Array.Copy(wholes, digitsVals, wholes.Length);
                return digitsVals;
            }
            else
                return wholes;
        }
        /// <summary>
        /// Gets the digis that make up the whole part of the number
        /// </summary>
        /// <param name="val">number to get the digits</param>
        /// <returns>array of digits form right to left from the decimal point</returns>
        public static int[] WholeDigits(float val)
        {
            return WholeDigits((double)val);
        }
        /// <summary>
        /// Gets the digis that make up the fraction part of the number
        /// </summary>
        /// <param name="val">number to get the digits</param>
        /// <returns>array of digits form left to right for the decimal point</returns>
        public static int[] FractionDigits(float val)
        {
            return FractionDigits((double)val);
        }
        /// <summary>
        /// Gets the digis that make up the whole part of the number
        /// </summary>
        /// <param name="val">number to get the digits</param>
        /// <returns>array of digits form right to left from the decimal point</returns>
        public static int[] WholeDigits(double val)
        {
            List<int> digits = new List<int>();
            int Whole = (int)val;
            while (Whole > 0)
            {
                digits.Add(Whole % 10);
                Whole /= 10;
            }
            return digits.ToArray();
        }
        /// <summary>
        /// Gets the digis that make up the fraction part of the number
        /// </summary>
        /// <param name="val">number to get the digits</param>
        /// <returns>array of digits form left to right for the decimal point</returns>
        public static int[] FractionDigits(double val)
        {
            List<int> digits = new List<int>();
            int Whole = (int)val;
            double Fraction = val - (double)Whole;
            Fraction *= 10;
            Whole = (int)Fraction % 10;
            while (Whole > 0)
            {
                digits.Add(Whole);
                Fraction *= 10;
                Whole = (int)Fraction % 10;
            }
            return digits.ToArray();
        }
    }
    /// <summary>
    /// Calculates a 8-bit fletcher checksum
    /// </summary>
    public class FletcherChecksum
    {
        public static void Calc(byte[] Packet, ref byte CS1, ref byte CS2)
        {
            CS1 = 0;
            CS2 = 0;
            for (int i = 0; i < Packet.Length; i++)
            {
                CS1 = (byte)(CS1 + Packet[i]);
                CS2 = (byte)(CS2 + CS1);
            }
        }
    }

    /// <summary>
    /// Calculates a 16-bit checksum, The checksum is an addition of all the data
    /// </summary>
    public class Checksum16Bit
    {
        public static ushort Calc(byte[] Packet)
        {
            ushort chksum = 0;
            for (int i = 0; i < Packet.Length; i++)
            {
                chksum += Packet[i];
            }
            return chksum;
        }
    }

    /// <summary>
    /// Calculates a 8-bit checksum. The checksum is an XOR of the data
    /// </summary>
    public class Checksum8Bit
    {
        /// <summary>
        /// String type to calculate the checksum on
        /// </summary>
        public enum StringType
        {
            ansi,
            unicode,
        }
        /// <summary>
        /// Calculate the checksum using the byte array
        /// </summary>
        /// <param name="Packet">array to checksum</param>
        /// <returns>checksum value</returns>
        public static byte Calc(byte[] Packet)
        {
            byte chksum = 0;
            for (int i = 0; i < Packet.Length; i++)
            {
                chksum ^= Packet[i];
            }
            return chksum;
        }
        /// <summary>
        /// Calculate the checksum using the string based on the string type
        /// </summary>
        /// <param name="Packet">string to checksum</param>
        /// <param name="type">string type</param>
        /// <returns>checksum value</returns>
        public static byte Calc(string Packet, StringType type)
        {
            byte chksum = 0;
            if (type == StringType.ansi)
            {
                for (int i = 0; i < Packet.Length; i++)
                {
                    chksum ^= (byte)Packet[i];
                }
                return chksum;
            }
            else
            {
                return Calc(Converter.GetBytes(Packet));
            }
        }
    }
    /// <summary>
    /// Static methods used for logging data
    /// </summary>
    public class Logging
    {
        /// <summary>
        /// Log to the application event log
        /// </summary>
        /// <param name="ErrorDescription">text to log</param>
        public static void LogException(string ErrorDescription)
        {
            // The name of our log in the event logs
            string Log = "Application";

            // Check to see fi the log for AspNetError exists on the machine
            //          If not, create it
            if ((!(EventLog.SourceExists(Log))))
            {
                EventLog.CreateEventSource(Log, Log);
            }

            // Now insert your exception information into the AspNetError event log
            EventLog logEntry = new EventLog();
            logEntry.Source = Log;
            logEntry.WriteEntry(ErrorDescription, EventLogEntryType.Error);
        }
    }

    /// <summary>
    /// Used to do a deep copy of an object
    /// http://www.codeproject.com/KB/cs/Deep_copy_of_objects.aspx
    /// </summary>
    public class Copy
    {
        // deep copy in separeate memory space
        /// <summary>
        /// Deep copy of an object
        /// </summary>
        /// <typeparam name="T">Type of the object</typeparam>
        /// <param name="obj">The object to copy</param>
        /// <returns>copied object</returns>
        public static T DeepCopy<T>(T obj) {
            if(obj==null)
                throw new ArgumentNullException("Object cannot be null");
            return (T)Process(obj);
        }

        static object Process(object obj) {
            if(obj==null)
                return null;
            Type type=obj.GetType();
            if(type.IsValueType || type==typeof(string)) {
                return obj;
            }
            else if(type.IsArray) {
                Type elementType=Type.GetType(
                     type.FullName.Replace("[]",string.Empty));
                var array=obj as Array;
                Array copied=Array.CreateInstance(elementType,array.Length);
                for(int i=0; i<array.Length; i++) {
                    copied.SetValue(Process(array.GetValue(i)),i);
                }
                return Convert.ChangeType(copied,obj.GetType());
            }
            else if(type.IsClass) {
                object toret=Activator.CreateInstance(obj.GetType());
                FieldInfo[] fields=type.GetFields(BindingFlags.Public| 
                            BindingFlags.NonPublic|BindingFlags.Instance);
                foreach(FieldInfo field in fields) {
                    object fieldValue=field.GetValue(obj);
                    if(fieldValue==null)
                        continue;
                    field.SetValue(toret,Process(fieldValue));
                }
                return toret;
            }
            else
                throw new ArgumentException("Unknown type");
        }
    }
    /// <summary>
    /// Static methods used to convert data in different endian formats
    /// </summary>
    public class Endian
    {
        /// <summary>
        /// The endian type
        /// </summary>
        public enum EndianType
        {
            Little,
            Big
        }

        /// <summary>
        /// Fix endianness of value
        /// </summary>
        /// <param name="val">the value to fix</param>
        /// <param name="EndianT">endianness to fix the value to</param>
        /// <returns></returns>
        public static double Fix(double val, EndianType EndianT)
        {
            if (((EndianT == EndianType.Little) && BitConverter.IsLittleEndian) || 
                ((EndianT == EndianType.Big) && !BitConverter.IsLittleEndian))
            {
                return val;
            }
            else
            {
                byte[] ValBytes = BitConverter.GetBytes(val);
                Array.Reverse(ValBytes);
                return BitConverter.ToDouble(ValBytes, 0);
            }
        }

        /// <summary>
        /// Fix endianness of value
        /// </summary>
        /// <param name="val">the value to fix</param>
        /// <param name="EndianT">endianness to fix the value to</param>
        /// <returns></returns>
        public static float Fix(float val, EndianType EndianT)
        {
            if (((EndianT == EndianType.Little) && BitConverter.IsLittleEndian) ||
                ((EndianT == EndianType.Big) && !BitConverter.IsLittleEndian))
            {
                return val;
            }
            else
            {
                byte[] ValBytes = BitConverter.GetBytes(val);
                Array.Reverse(ValBytes);
                return BitConverter.ToSingle(ValBytes, 0);
            }
        }

        /// <summary>
        /// Fix endianness of value
        /// </summary>
        /// <param name="val">the value to fix</param>
        /// <param name="EndianT">endianness to fix the value to</param>
        /// <returns></returns>
        public static short Fix(short val, EndianType EndianT)
        {
            if (((EndianT == EndianType.Little) && BitConverter.IsLittleEndian) ||
                ((EndianT == EndianType.Big) && !BitConverter.IsLittleEndian))
            {
                return val;
            }
            else
            {
                byte[] ValBytes = BitConverter.GetBytes(val);
                Array.Reverse(ValBytes);
                return BitConverter.ToInt16(ValBytes, 0);
            }
        }

        /// <summary>
        /// Fix endianness of value
        /// </summary>
        /// <param name="val">the value to fix</param>
        /// <param name="EndianT">endianness to fix the value to</param>
        /// <returns></returns>
        public static ushort Fix(ushort val, EndianType EndianT)
        {
            if (((EndianT == EndianType.Little) && BitConverter.IsLittleEndian) ||
                ((EndianT == EndianType.Big) && !BitConverter.IsLittleEndian))
            {
                return val;
            }
            else
            {
                byte[] ValBytes = BitConverter.GetBytes(val);
                Array.Reverse(ValBytes);
                return BitConverter.ToUInt16(ValBytes, 0);
            }
        }

        /// <summary>
        /// Fix endianness of value
        /// </summary>
        /// <param name="val">the value to fix</param>
        /// <param name="EndianT">endianness to fix the value to</param>
        /// <returns></returns>
        public static int Fix(int val, EndianType EndianT)
        {
            if (((EndianT == EndianType.Little) && BitConverter.IsLittleEndian) ||
                ((EndianT == EndianType.Big) && !BitConverter.IsLittleEndian))
            {
                return val;
            }
            else
            {
                byte[] ValBytes = BitConverter.GetBytes(val);
                Array.Reverse(ValBytes);
                return BitConverter.ToInt32(ValBytes, 0);
            }
        }

        /// <summary>
        /// Fix endianness of value
        /// </summary>
        /// <param name="val">the value to fix</param>
        /// <param name="EndianT">endianness to fix the value to</param>
        /// <returns></returns>
        public static uint Fix(uint val, EndianType EndianT)
        {
            if (((EndianT == EndianType.Little) && BitConverter.IsLittleEndian) ||
                ((EndianT == EndianType.Big) && !BitConverter.IsLittleEndian))
            {
                return val;
            }
            else
            {
                byte[] ValBytes = BitConverter.GetBytes(val);
                Array.Reverse(ValBytes);
                return BitConverter.ToUInt32(ValBytes, 0);
            }
        }
    }

    /// <summary>
    /// Converts structs to byte arrays and vice versa
    /// http://bytes.com/groups/net-c/236808-how-convert-structure-byte-array
    /// </summary>
    public class Converter
    {
        /// <summary>
        /// Converts a structure to a byte array
        /// </summary>
        /// <param name="obj">struct object</param>
        /// <returns>array of bytes</returns>
        public static byte[] StructureToByteArray(object obj)
        {
            int len = Marshal.SizeOf(obj);
            byte[] arr = new byte[len];
            IntPtr ptr = Marshal.AllocHGlobal(len);
            Marshal.StructureToPtr(obj, ptr, true);
            Marshal.Copy(ptr, arr, 0, len);
            Marshal.FreeHGlobal(ptr);
            return arr;
        }

        /// <summary>
        /// Converts a byte array to a structure
        /// </summary>
        /// <param name="bytearray">byte array of structure</param>
        /// <param name="obj">ref object of struct</param>
        public static void ByteArrayToStructure(byte[] bytearray, ref object obj)
        {
            int len = Marshal.SizeOf(obj);
            IntPtr i = Marshal.AllocHGlobal(len);
            Marshal.Copy(bytearray, 0, i, len);
            obj = Marshal.PtrToStructure(i, obj.GetType());
            Marshal.FreeHGlobal(i);
        }

        /// <summary>
        /// Converts unicode string to byte array. First 4 bytes holds the int length of the string
        /// </summary>
        /// <param name="text">Unicode string to convert</param>
        /// <returns>byte array of string</returns>
        public static byte[] GetBytes(string text)
        {
            List<byte> StrPck = new List<byte>();
            StrPck.AddRange(BitConverter.GetBytes(text.Length));
            StrPck.AddRange(ASCIIEncoding.Unicode.GetBytes(text));
            return StrPck.ToArray();
        }

        /// <summary>
        /// Converts a unicode byte array to a string. First 4 bytes holds the int length of the string
        /// </summary>
        /// <param name="bytes">byte array of the string</param>
        /// <returns>string</returns>
        public static string GetString(byte[] bytes)
        {
            int len = BitConverter.ToInt32(bytes, 0);
            return new System.Text.UnicodeEncoding().GetString(bytes, 4, len*2);
        }
    }

    /// <summary>
    /// Class used for execution performance calculation
    /// </summary>
    public class Performance
    {
        Stopwatch ExecutionTime = new Stopwatch();

        /// <summary>
        /// Starts the timers
        /// </summary>
        public void Start()
        {
            ExecutionTime.Reset();
            ExecutionTime.Start();
        }

        /// <summary>
        /// Stop the timers
        /// </summary>
        /// <param name="text">string to output on the debug window</param>
        public void Stop(string text)
        {
            ExecutionTime.Stop();
#if DEBUG
            Debug.WriteLine(text);
            Debug.Indent();
            Debug.WriteLine(ExecutionTime.Elapsed.TotalMilliseconds.ToString() + " ms");
            Debug.Unindent();
#endif
        }
    }
    /// <summary>
    /// Static methods to assist with trigometric calculations
    /// </summary>
    public class TrigCalc
    { 
        /// <summary>
        /// Calculates the length of the side opposite the reference angle.
        /// </summary>
        /// <param name="adjacent">the length of the side adjacent to the reference angle</param>
        /// <param name="angle">the size of the reference angle</param>
        /// <returns>the length of the side opposite the reference angle.</returns>
        public static double OppositeTan(double adjacent, double angle)
        { 
            return adjacent * Math.Tan(math.Deg2Rad(angle));
        }


        /// <summary>
        /// Calculates the length of the side adjacent to the reference angle.
        /// </summary>
        /// <param name="hyp">the length on the hypotenuse</param>
        /// <param name="angle">the size of the reference angle</param>
        /// <returns>the length of the side adjacent to the reference angle.</returns>
        public static double AdjacentCos(double hyp, double angle)
        {
            return hyp * Math.Cos(math.Deg2Rad(angle));
        }

        /// <summary>
        /// Calculates the length of the side opposite to the reference angle.
        /// </summary>
        /// <param name="hyp">the length on the hypotenuse</param>
        /// <param name="angle">the size of the reference angle</param>
        /// <returns>the length of the side opposite the reference angle.</returns>
        public static double OppositeSin(double hyp, double angle)
        {
            return hyp * Math.Sin(math.Deg2Rad(angle));
        }
    }
    /// <summary>
    /// Static methods for optics calculations
    /// </summary>
    public class Optics
    {
        /// <summary>
        /// Used to determine the image sensor width given the focal length and the field of view
        /// </summary>
        /// <param name="FocalLength">focal length value</param>
        /// <param name="FieldOfView">field of view value</param>
        /// <returns>width</returns>
        public static double ImageSensorWidth(double FocalLength, double FieldOfView)
        {
            double Width = 2 * Math.Tan(math.Deg2Rad(FieldOfView / 2)) * FocalLength;
            return Width;
        }

        /// <summary>
        /// Used to determine the field of view in degrees given the sensor width and the focal length
        /// </summary>
        /// <param name="SensorWidth">the image sensor width</param>
        /// <param name="FocalLength">the focal length value</param>
        /// <returns>the field of view in degrees</returns>
        public static double FieldOfView(double SensorWidth, double FocalLength)
        {
            double FOV = 2 * Math.Atan(SensorWidth / (2 * FocalLength));
            return math.Rad2Deg(FOV);
        }
    }
    /// <summary>
    /// Static methods for temperature conversion
    /// </summary>
    public class Temperature
    {
        /// <summary>
        /// Degrees Celcius to Fahrenheit
        /// </summary>
        /// <param name="temperature">The temperature in degrees</param>
        /// <returns>The temperature in fahrenheit</returns>
        public static double CelsiusToFahrenheit(double temperature)
        {
            return (5/9) * (temperature-32);
        }
        /// <summary>
        /// Degrees Fahrenheit to Celcius
        /// </summary>
        /// <param name="temperature">The temperature in fahrenheit</param>
        /// <returns>The temperature in degrees</returns>
        public static double FahrenheitToCelsius(double temperature)
        {
            return (9/5)*temperature+32;
        }

    }
    /// <summary>
    /// Static methods for distance calculations
    /// </summary>
    public class Distance
    {
        /// <summary>
        /// Convert kilometers to miles
        /// </summary>
        /// <param name="kilometers">distance in kilometers to convert</param>
        /// <returns>distance calculated in miles</returns>
        public static double KmToMiles(double kilometers)
        {
            return kilometers * 0.621371192;
        }
        /// <summary>
        /// Convert miles to kilometers
        /// </summary>
        /// <param name="miles">distance in miles to convert</param>
        /// <returns>distance calculated in kilometers</returns>
        public static double MilesToKm(double miles)
        {
            return miles * 1.609344;
        }
        /// <summary>
        /// Convert meters to feet
        /// </summary>
        /// <param name="meters">distance in meters to convert</param>
        /// <returns>distance calculated in feet</returns>
        public static double MetersToFeet(double meters)
        {
            return meters * 3.2808;
        }
        /// <summary>
        /// Convert feet to meters
        /// </summary>
        /// <param name="feet">distance in feet to convert</param>
        /// <returns>distance calculated in meters</returns>
        public static double FeetToMeters(double feet)
        {
            return feet * 0.3048;
        }
        /// <summary>
        /// Convert nautical miles to meters
        /// </summary>
        /// <param name="NM">distance in nautical miles to convert</param>
        /// <returns>distance calculated in meters</returns>
        public static double NauticalMilesToMeters(double NM)
        {
            return NM * 1852;
        }
        /// <summary>
        /// Convert meters to nautical miles
        /// </summary>
        /// <param name="meters">distance in meters to convert</param>
        /// <returns>distance calculated in nautical miles</returns>
        public static double MetersToNauticalMiles(double meters)
        {
            return meters / 1852;
        }
    }

    //public class DataBase
    //{
    //    public static DataTable BindToEnum(Type enumType)
    //    {
    //        string[] names = Enum.GetNames(enumType);
    //        int[] values = (int[])Enum.GetValues(enumType);
    //        DataTable dt = new DataTable();
    //        dt.Columns.Add("Key", typeof(string));
    //        dt.Columns.Add("Value", typeof(int));
    //        for (int i = 0; i < names.Length; i++)
    //        {
    //            //DataRow dr = new 
    //            dr["Key"] = names[i];
    //            dr["Value"] = values[i];
    //            dt.Rows.Add(dr);
    //        }
    //        return dt;
    //    }
    //}
}
