using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.IO.Ports;
using System.IO;
using System.Diagnostics;
//using ZedGraph;
using System.Runtime.InteropServices;

namespace WindowsFormsApplication1
{
    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    internal struct IMU_Data_Struct
    {
        public float PLACEHOLDER; // This needs to go here
        public UInt16 gyroX1;
        public UInt16 gyroY1;
        public UInt16 gyroZ1;
        public UInt16 magX1;
        public UInt16 magY1;
        public UInt16 magZ1;
        public UInt16 accX1;
        public UInt16 accY1;
        public UInt16 accZ1;
        
        public byte IMUtemp1;
        public char adcTemp1;

        public UInt16 gyroX2;
        public UInt16 gyroY2;
        public UInt16 gyroZ2;
        public UInt16 magX2;
        public UInt16 magY2;
        public UInt16 magZ2;
        public UInt16 accX2;
        public UInt16 accY2;
        public UInt16 accZ2;

        public byte IMUtemp2;
        public char adcTemp2;

        public UInt16 gyroX3;
        public UInt16 gyroY3;
        public UInt16 gyroZ3;
        public UInt16 magX3;
        public UInt16 magY3;
        public UInt16 magZ3;
        public UInt16 accX3;
        public UInt16 accY3;
        public UInt16 accZ3;

        public byte IMUtemp3;
        public char adcTemp3;

        public UInt16 gyroX4;
        public UInt16 gyroY4;
        public UInt16 gyroZ4;
        public UInt16 magX4;
        public UInt16 magY4;
        public UInt16 magZ4;
        public UInt16 accX4;
        public UInt16 accY4;
        public UInt16 accZ4;

        public byte IMUtemp4;
        public char adcTemp4;

        public byte success;
       
    };

    public partial class Form1 : Form
    {
        //************************************* VARIABLE DECLARATIONS *************************************
        HDLC txDataFrame;
        HDLC rxDataFrame;
        int baudR = 230400;         // Baud rate
        string RXString;
        string pc_logfile_name = "";
        bool pc_logfile_open = false;
        bool log_data_checked_old = false;

        int count = 0;
        int ValUpdate = 60; 

        Stopwatch sw = new Stopwatch();
        //StreamWriter log; // For the logfile
        BinaryWriter log;
        IMU_Data_Struct IMU_data;

        // Constructor for the form
        public Form1()
        {
            txDataFrame = new HDLC();
            rxDataFrame = new HDLC();
            InitializeComponent();

            // Add the KEYDOWN event handler to handle pressed keys
            KeyDown += new KeyEventHandler(Form1_KeyDown);
            this.KeyPreview = true;

            // Add the FORM CLOSE event handler
            this.Closing += new CancelEventHandler(Form_Closing);
        }

        //************************************* EVENT HANDLERS *************************************

        // Event handler for KEYDOWN action
        void Form1_KeyDown(object sender, KeyEventArgs e)
        {
            switch (e.KeyCode)
            {
                case Keys.A: //steerLeft();
                    break;
                case Keys.D: //steerRight();
                    break;
                case Keys.W: //increaseThrottle();
                    break;
                case Keys.S: //decreaseThrottle();
                    break;
                case Keys.Space: //stopCar();
                    break;
                case Keys.Q: //centerSteering();
                    break;
                case Keys.K: //killCar();
                    break;
                case Keys.R: //reviveCar();
                    break;
                case Keys.P: //togglePing();
                    break;
                case Keys.M: //switchToManual();
                    break;
                case Keys.N: //swtichToAuto();
                    break;
            }
        }

        // Event handler for when the SERIAL PORT OPEN BUTTON is clicked
        private void openButton_Click(object sender, EventArgs e)
        {
            try
            {
                if (!(serialPort1.IsOpen))
                {
                    // Start the STOPWATCH
                    sw.Start();
                    //time_old = sw.ElapsedMilliseconds;

                    // Get the Port number from the Combo Box
                    serialPort1.PortName = serialComboBox.Text;

                    // Get the Baud Rate from the other Combo Box
                    serialPort1.BaudRate = baudR;

                    // Open the Serial Port
                    serialPort1.Open();

                    // Change the text
                    openButton.Text = "Close Port";
                    serialComboBox.Enabled = false;

                    serialPort1.DataReceived += new SerialDataReceivedEventHandler(serialDataReceivedHandler);
                }

                else if (serialPort1.IsOpen)
                {
                    
                    // Close serial port
                    serialPort1.Close();
                    openButton.Text = "Disabled";
                    // Disable button
                    openButton.Enabled = false;
                    // Close the filestream if it is open
                    if (pc_logfile_open == true)
                        log.Close();
                }
            }

            catch (Exception ex)
            {
                MessageBox.Show(ex.Message);
                MessageBox.Show("2");
            }
        }

        private void serialDataReceivedHandler(object sender, SerialDataReceivedEventArgs e)
        {
            byte[] originalFrame = new byte[512];//this was 512
           
            while (serialPort1.IsOpen && (serialPort1.BytesToRead > 0))
            {
                try
                {
                    int j = 0;
                    byte[] Data = rxDataFrame.ProcessRxByte((byte)serialPort1.ReadByte(), originalFrame);
                    int P = 0;
                    if (Data != null)
                    {
                        //DATA LOG Data EVERY TIME
                        if (pc_logfile_open == true)
                        {
                            byte[]  result = new byte[200];
                            Array.Copy(originalFrame, 0, result, 0, 200);
                            log.Write(result);
                            //log.Write("\r\n");
                        }
                        //RUN THIS EVERY X TIMES
                        count = count + 1;
                        if (count == ValUpdate)
                        {
                            UInt16 command = System.BitConverter.ToUInt16(Data, 0);
                            if (command == 0x0101)
                            {
                                IMU_data = ByteCaster.GetStructure<IMU_Data_Struct>(Data);
                                this.Invoke(new EventHandler(updateGUI_Information));
                            }
                            count = 0;
                        }
                    }
                    else
                    {
                        this.Invoke(new EventHandler(updateGUI_Information));
                        
                    }
                }
                catch (Exception h)
                {
                    // Close the filestream if it is open
                    if (pc_logfile_open == true)
                        log.Close();
                    MessageBox.Show(h.Message);
                }
            }
        }

        
        private void logMethod()
        {
            try
            {
               
                if (log_data_checked_old != checkBoxLogging.Checked)
                {
                    log_data_checked_old = checkBoxLogging.Checked;

                    // Enable the logging
                    if (checkBoxLogging.Checked == true)
                    {
                        // Logfile filename is: DIMA_'current date'_'current time'.txt
                        pc_logfile_name = "/Users/Callen Fisher/Desktop/logdata/" + getFormattedDateTime() + ".dat";

                        // Create the PC logfile
                        if (!File.Exists(pc_logfile_name))
                        {
                            log = new BinaryWriter(File.Open(pc_logfile_name, FileMode.Create));
                            //log = new StreamWriter(pc_logfile_name);
                            pc_logfile_open = true;
                        }
                        // ADD ANY COLUMNS TO LOG
                        //log.WriteLine("GyroX1,GyroY1,GyroZ1,AccX1,AccY1,AccZ1,MagX1,MagY1,MagZ1,Temp1,GyroX2,GyroY2,GyroZ2,AccX2,AccY2,AccZ2,MagX2,MagY2,MagZ2,Temp2,GyroX3,GyroY3,GyroZ3,AccX3,AccY3,AccZ3,MagX3,MagY3,MagZ3,Temp3,GyroX4,GyroY4,GyroZ4,AccX4,AccY4,AccZ4,MagX4,MagY4,MagZ4,Temp4");

                        // First byte is NON-ZERO
                        // 12 characters containing the filename. Month/Day/Hours/Minutes/ - "08091541.dat"

                        int month = DateTime.Now.Month;
                        int day = DateTime.Now.Day;
                        int hours = DateTime.Now.Hour;
                        int minutes = DateTime.Now.Minute;
                        //int seconds = DateTime.Now.Second;

                        StringBuilder sb = new StringBuilder(19);

                        // Add the MONTH
                        if (month < 10)
                            sb.Append("0");
                        sb.Append(month.ToString());

                        // Add the DAY
                        if (day < 10)
                            sb.Append("0");

                        sb.Append(day.ToString());

                        // Add the HOURS
                        if (hours < 10)
                            sb.Append("0");

                        sb.Append(hours.ToString());

                        // Add the MINUTES
                        if (minutes < 10)
                            sb.Append("0");

                        sb.Append(minutes.ToString());                     
                    }
                    // Disable the logging 
                    else
                    {
                        // Close the PC logfile
                        pc_logfile_open = false;
                        log.Close();
                        MessageBox.Show("file closed");
                    }
                }
                
            }

            catch (Exception ex)
            {
                // Close the filestream if it is open
                if (pc_logfile_open == true)
                    log.Close();

                MessageBox.Show(ex.Message);
                MessageBox.Show("program crashed - please restart");
            }
        }
        // Update telemetry data in GUI as well as write the data to the logfile
        private void updateGUI_Information(object sender, EventArgs e)
        {
            if (count == ValUpdate)
            {

                //TODO twos comp and scaling factor

                float accScale = 0.000786F;
                float magScale = 0.0043F;
                float magScale2 = 0.0049F;
                float gyroScale = 0.07F;
		float tempScale=1;

                String gyroX1 = (twosComp(IMU_data.gyroX1)*gyroScale).ToString();
                String gyroY1 = (twosComp(IMU_data.gyroY1)*gyroScale).ToString();
                String gyroZ1 = (twosComp(IMU_data.gyroZ1)*gyroScale).ToString();
                String accX1 = (twosComp(IMU_data.accX1)*accScale).ToString();
                String accY1 = (twosComp(IMU_data.accY1)*accScale).ToString();
                String accZ1 = (twosComp(IMU_data.accZ1)*accScale).ToString();
                String magX1 = (twosComp(IMU_data.magX1)*magScale).ToString();
                String magY1 = (twosComp(IMU_data.magY1)*magScale).ToString();
                String magZ1 = (twosComp(IMU_data.magZ1)*magScale2).ToString();
                String IMUTemp1 = (twosCompTemp(IMU_data.IMUtemp1)*tempScale).ToString();
                String ADCtemp1 = IMU_data.adcTemp1.ToString();

                String gyroX2 = (twosComp(IMU_data.gyroX2)*gyroScale).ToString();
                String gyroY2 = (twosComp(IMU_data.gyroY2)*gyroScale).ToString();
                String gyroZ2 = (twosComp(IMU_data.gyroZ2)*gyroScale).ToString();
                String accX2 = (twosComp(IMU_data.accX2)*accScale).ToString();
                String accY2 = (twosComp(IMU_data.accY2)*accScale).ToString();
                String accZ2 = (twosComp(IMU_data.accZ2)*accScale).ToString();
                String magX2 = (twosComp(IMU_data.magX2)*magScale).ToString();
                String magY2 = (twosComp(IMU_data.magY2)*magScale).ToString();
                String magZ2 = (twosComp(IMU_data.magZ2)*magScale2).ToString();
                String IMUTemp2 = (twosCompTemp(IMU_data.IMUtemp2)*tempScale).ToString();//TODO 
                String ADCtemp2 = IMU_data.adcTemp2.ToString();

                String gyroX3 = (twosComp(IMU_data.gyroX3)*gyroScale).ToString();
                String gyroY3 = (twosComp(IMU_data.gyroY3)*gyroScale).ToString();
                String gyroZ3 = (twosComp(IMU_data.gyroZ3)*gyroScale).ToString();
                String accX3 = (twosComp(IMU_data.accX3)*accScale).ToString();
                String accY3 = (twosComp(IMU_data.accY3)*accScale).ToString();
                String accZ3 = (twosComp(IMU_data.accZ3)*accScale).ToString();
                String magX3 = (twosComp(IMU_data.magX3)*magScale).ToString();
                String magY3 = (twosComp(IMU_data.magY3)*magScale).ToString();
                String magZ3 = (twosComp(IMU_data.magZ3)*magScale2).ToString();
                String IMUTemp3 = (twosCompTemp(IMU_data.IMUtemp3)*tempScale).ToString();
                String ADCtemp3 = IMU_data.adcTemp3.ToString();

                String gyroX4 = (twosComp(IMU_data.gyroX4)*gyroScale).ToString();
                String gyroY4 = (twosComp(IMU_data.gyroY4)*gyroScale).ToString();
                String gyroZ4 = (twosComp(IMU_data.gyroZ4)*gyroScale).ToString();
                String accX4 = (twosComp(IMU_data.accX4)*accScale).ToString();
                String accY4 = (twosComp(IMU_data.accY4)*accScale).ToString();
                String accZ4 = (twosComp(IMU_data.accZ4)*accScale).ToString();
                String magX4 = (twosComp(IMU_data.magX4)*magScale).ToString();
                String magY4 = (twosComp(IMU_data.magY4)*magScale).ToString();
                String magZ4 = (twosComp(IMU_data.magZ4)*magScale2).ToString();
                String IMUTemp4 = (twosCompTemp(IMU_data.IMUtemp4)*tempScale).ToString();
                String ADCtemp4 = IMU_data.adcTemp4.ToString();

                // Update the information in the GUI
                labelGyroXVal.Text = gyroX1;
                labelGyroYVal.Text = gyroY1;
                labelGyroZVal.Text = gyroZ1;
                labelAccXVal.Text = accX1;
                labelAccYVal.Text = accY1;
                labelAccZVal.Text = accZ1;
                labelMagXVal.Text = magX1;
                labelMagYVal.Text = magY1;
                labelMagZVal.Text = magZ1;
                labelIMUTempVal.Text = IMUTemp1;
                ADCtempLABEL1.Text = ADCtemp1;

                labelGyroXVal2.Text = gyroX2;
                labelGyroYVal2.Text = gyroY2;
                labelGyroZVal2.Text = gyroZ2;
                labelAccXVal2.Text = accX2;
                labelAccYVal2.Text = accY2;
                labelAccZVal2.Text = accZ2;
                labelMagXVal2.Text = magX2;
                labelMagYVal2.Text = magY2;
                labelMagZVal2.Text = magZ2;
                labelIMUTempVal2.Text = IMUTemp2;
                ADCtempLABEL2.Text = ADCtemp2;

                labelGyroXVal3.Text = gyroX3;
                labelGyroYVal3.Text = gyroY3;
                labelGyroZVal3.Text = gyroZ3;
                labelAccXVal3.Text = accX3;
                labelAccYVal3.Text = accY3;
                labelAccZVal3.Text = accZ3;
                labelMagXVal3.Text = magX3;
                labelMagYVal3.Text = magY3;
                labelMagZVal3.Text = magZ3;
                labelIMUTempVal3.Text = IMUTemp3;
                ADCtempLABEL3.Text = ADCtemp3;

                labelGyroXVal4.Text = gyroX4;
                labelGyroYVal4.Text = gyroY4;
                labelGyroZVal4.Text = gyroZ4;
                labelAccXVal4.Text = accX4;
                labelAccYVal4.Text = accY4;
                labelAccZVal4.Text = accZ4;
                labelMagXVal4.Text = magX4;
                labelMagYVal4.Text = magY4;
                labelMagZVal4.Text = magZ4;
                labelIMUTempVal4.Text = IMUTemp4;
                ADCtempLABEL4.Text = ADCtemp4;

                // Update CRC and ESC CHAR counters
                labelCRCFailsVal.Text = rxDataFrame.getNumberOfCRCFails().ToString();
                labelEscCharFailsVal.Text = rxDataFrame.getNumberOfescapeCharFails().ToString();


                successByteLabel.Text = IMU_data.success.ToString();
                //// Construct a string to be printed to the logfile
                //StringBuilder logString = new StringBuilder(1000);
                //logString.Append(gyroX1 + ";");
                //logString.Append(gyroY1 + ";");
                //logString.Append(gyroZ1 + ";");
                //logString.Append(accX1 + ";");
                //logString.Append(accY1 + ";");
                //logString.Append(accZ1 + ";");
                //logString.Append(magX1 + ";");
                //logString.Append(magY1 + ";");
                //logString.Append(magZ1 + ";");
                //logString.Append(IMUTemp1 + ";");

                //logString.Append(gyroX2 + ";");
                //logString.Append(gyroY2 + ";");
                //logString.Append(gyroZ2 + ";");
                //logString.Append(accX2 + ";");
                //logString.Append(accY2 + ";");
                //logString.Append(accZ2 + ";");
                //logString.Append(magX2 + ";");
                //logString.Append(magY2 + ";");
                //logString.Append(magZ2 + ";");
                //logString.Append(IMUTemp2 + ";");

                //logString.Append(gyroX3 + ";");
                //logString.Append(gyroY3 + ";");
                //logString.Append(gyroZ3 + ";");
                //logString.Append(accX3 + ";");
                //logString.Append(accY3 + ";");
                //logString.Append(accZ3 + ";");
                //logString.Append(magX3 + ";");
                //logString.Append(magY3 + ";");
                //logString.Append(magZ3 + ";");
                //logString.Append(IMUTemp3 + ";");

                //logString.Append(gyroX4 + ";");
                //logString.Append(gyroY4 + ";");
                //logString.Append(gyroZ4 + ";");
                //logString.Append(accX4 + ";");
                //logString.Append(accY4 + ";");
                //logString.Append(accZ4 + ";");
                //logString.Append(magX4 + ";");
                //logString.Append(magY4 + ";");
                //logString.Append(magZ4 + ";");
                //logString.Append(IMUTemp4);
                //count = 0;
                //if (pc_logfile_open == true)
                //    log.WriteLine(logString.ToString());
            }
        }
        float twosComp(UInt16 v)
        {
            //u16 v=val;
	        float temp=0;
            UInt16 one=5;
            UInt16 two=5;
            UInt16  tempTest=(UInt16)(one & two);

	        if( (UInt16)(v & (UInt16)(32768))==32768 )
	        {
		        temp=-32768;
	        }
	        if((UInt16)(v & 1)==1)
	        {
		        temp=temp+1;
	        }
	        if((UInt16)(v & 2)==2)
	        {
		        temp=temp+2;
	        }
	        if((UInt16)(v & 4)==4)
	        {
		        temp=temp+4;
	        }
	        if((UInt16)(v & 8)==8)
	        {
		        temp=temp+8;
	        }
	        if((UInt16)(v & 16)==16)
	        {
		        temp=temp+16;
	        }
	        if((UInt16)(v & 32)==32)
	        {
		        temp=temp+32;
	        }
	        if((UInt16)(v & 64)==64)
	        {
		        temp=temp+64;
	        }
	        if((UInt16)(v & 128)==128)
	        {
		        temp=temp+128;
	        }
	        if((UInt16)(v & 256)==256)
	        {
		        temp=temp+256;
	        }
	        if((UInt16)(v & 512)==512)
	        {
		        temp=temp+512;
	        }
	        if((UInt16)(v & 1024)==1024)
	        {
		        temp=temp+1024;
	        }
	        if((UInt16)(v & 2048)==2048)
	        {
	        	temp=temp+2048;
	        }
	        if((UInt16)(v & 4096)==4096)
	        {
		        temp=temp+4096;
	        }
	        if((UInt16)(v & 8192)==8192)
	        {
		        temp=temp+8192;
	        }
	        if((UInt16)(v & 16384)==16384)
	        {
		        temp=temp+16384;
	        }

	    return temp;
        }
float twosCompTemp(UInt16 v)
        {
            float temp = 0;
            
            if ((UInt16)(v & (UInt16)(32768)) == 32768)
            {
                temp = -128;
            }
            if ((UInt16)(v & 1) == 1)
            {
                temp = (temp + 1);
            }
            if ((UInt16)(v & 2) == 2)
            {
                temp = (temp + 2);
            }
            if ((UInt16)(v & 4) == 4)
            {
                temp = (temp + 4);
            }
            if ((UInt16)(v & 8) == 8)
            {
                temp = (temp + 8);
            }
            if ((UInt16)(v & 16) == 16)
            {
                temp = (temp + 16);
            }
            if ((UInt16)(v & 32) == 32)
            {
                temp = (temp + 32);
            }
            if ((UInt16)(v & 64) == 64)
            {
                temp =( temp + 64);
            }
   

            return temp;
        }
        // Display text in DEBUG textbox
        private void DisplayText(object sender, EventArgs e)
        {
            textBoxDebug.AppendText(RXString);
        }

        // When form is closed, make sure serial ports are closed
        private void Form1_FormClosing(object sender, FormClosingEventArgs e)
        {
            if (serialPort1.IsOpen) serialPort1.Close();
        }
    
        private void Form_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            this.Invoke(new EventHandler(closeFormSafely));
        }

        private void closeFormSafely(object sender, EventArgs e)
        {
            // Close serial port
            if (serialPort1.IsOpen)
            {
                serialPort1.Close();
            }
        }

        private String getFormattedDateTime()
        {
            int year = DateTime.Now.Year;
            int month = DateTime.Now.Month;
            int day = DateTime.Now.Day;
            int hours = DateTime.Now.Hour;
            int minutes = DateTime.Now.Minute;
            int seconds = DateTime.Now.Second;

            StringBuilder sb = new StringBuilder(19);

            // Add the YEAR
            sb.Append(year.ToString());
            sb.Append("_");

            // Add the MONTH
            if (month < 10)
                sb.Append("0");

            sb.Append(month.ToString());
            sb.Append("_");

            // Add the DAY
            if (day < 10)
                sb.Append("0");

            sb.Append(day.ToString());
            sb.Append("_");

            // Add the HOURS
            if (hours < 10)
                sb.Append("0");

            sb.Append(hours.ToString());
            sb.Append("_");

            // Add the MINUTES
            if (minutes < 10)
                sb.Append("0");

            sb.Append(minutes.ToString());
            sb.Append("_");

            // Add the SECONDS
            if (seconds < 10)
                sb.Append("0");

            sb.Append(seconds.ToString());

            String formattedDateTime = sb.ToString();

            return formattedDateTime;
        }
        private void buttonClearDebug_Click(object sender, EventArgs e)
        {
            textBoxDebug.Clear();
        }

        // Function to determine if a bit in a half-word is set or not
        // pos=0 is LSb, pos=15 is MSb
        private bool IsBitSet(UInt16 b, int pos)
        {
            return (b & (1 << pos)) != 0;
        }
        private void buttonResetCounters_Click(object sender, EventArgs e)
        {
            rxDataFrame.resetCounters();
        }

        private void checkBoxLogging_CheckedChanged(object sender, EventArgs e)
        {
            //send sync command
            byte[] dataOut = new byte[1];
            dataOut[0] = 0x7E;
            if (serialPort1.IsOpen)
            {
                serialPort1.Write(dataOut, 0,1);
            }
            logMethod();
        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }
    }
}