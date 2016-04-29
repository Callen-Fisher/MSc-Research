namespace WindowsFormsApplication1
{
    partial class Form1
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.components = new System.ComponentModel.Container();
            this.serialComboBox = new System.Windows.Forms.ComboBox();
            this.labelSerialPort = new System.Windows.Forms.Label();
            this.openButton = new System.Windows.Forms.Button();
            this.serialPort1 = new System.IO.Ports.SerialPort(this.components);
            this.textBoxDebug = new System.Windows.Forms.TextBox();
            this.labelDebug = new System.Windows.Forms.Label();
            this.buttonClearDebug = new System.Windows.Forms.Button();
            this.textBoxWarnings = new System.Windows.Forms.TextBox();
            this.groupBoxWarnings = new System.Windows.Forms.GroupBox();
            this.groupBoxCRC = new System.Windows.Forms.GroupBox();
            this.buttonResetCounters = new System.Windows.Forms.Button();
            this.labelEscCharFailsVal = new System.Windows.Forms.Label();
            this.labelCRCFailsVal = new System.Windows.Forms.Label();
            this.labelEscCharFails = new System.Windows.Forms.Label();
            this.labelCRCFails = new System.Windows.Forms.Label();
            this.checkBoxLogging = new System.Windows.Forms.CheckBox();
            this.labelGyroXVal = new System.Windows.Forms.Label();
            this.labelGyroYVal = new System.Windows.Forms.Label();
            this.labelGyroZVal = new System.Windows.Forms.Label();
            this.labelMagXVal = new System.Windows.Forms.Label();
            this.labelMagYVal = new System.Windows.Forms.Label();
            this.labelMagZVal = new System.Windows.Forms.Label();
            this.labelAccXVal = new System.Windows.Forms.Label();
            this.labelAccYVal = new System.Windows.Forms.Label();
            this.labelAccZVal = new System.Windows.Forms.Label();
            this.labelIMUTempVal = new System.Windows.Forms.Label();
            this.label1 = new System.Windows.Forms.Label();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.label4 = new System.Windows.Forms.Label();
            this.panel1 = new System.Windows.Forms.Panel();
            this.ADCtempLABEL1 = new System.Windows.Forms.Label();
            this.panel2 = new System.Windows.Forms.Panel();
            this.ADCtempLABEL2 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label7 = new System.Windows.Forms.Label();
            this.label8 = new System.Windows.Forms.Label();
            this.labelIMUTempVal2 = new System.Windows.Forms.Label();
            this.labelAccZVal2 = new System.Windows.Forms.Label();
            this.labelAccYVal2 = new System.Windows.Forms.Label();
            this.labelAccXVal2 = new System.Windows.Forms.Label();
            this.labelMagZVal2 = new System.Windows.Forms.Label();
            this.labelMagYVal2 = new System.Windows.Forms.Label();
            this.labelMagXVal2 = new System.Windows.Forms.Label();
            this.labelGyroZVal2 = new System.Windows.Forms.Label();
            this.labelGyroYVal2 = new System.Windows.Forms.Label();
            this.labelGyroXVal2 = new System.Windows.Forms.Label();
            this.panel3 = new System.Windows.Forms.Panel();
            this.ADCtempLABEL3 = new System.Windows.Forms.Label();
            this.label19 = new System.Windows.Forms.Label();
            this.label20 = new System.Windows.Forms.Label();
            this.label21 = new System.Windows.Forms.Label();
            this.label22 = new System.Windows.Forms.Label();
            this.labelIMUTempVal3 = new System.Windows.Forms.Label();
            this.labelAccZVal3 = new System.Windows.Forms.Label();
            this.labelAccYVal3 = new System.Windows.Forms.Label();
            this.labelAccXVal3 = new System.Windows.Forms.Label();
            this.labelMagZVal3 = new System.Windows.Forms.Label();
            this.labelMagYVal3 = new System.Windows.Forms.Label();
            this.labelMagXVal3 = new System.Windows.Forms.Label();
            this.labelGyroZVal3 = new System.Windows.Forms.Label();
            this.labelGyroYVal3 = new System.Windows.Forms.Label();
            this.labelGyroXVal3 = new System.Windows.Forms.Label();
            this.panel4 = new System.Windows.Forms.Panel();
            this.ADCtempLABEL4 = new System.Windows.Forms.Label();
            this.label33 = new System.Windows.Forms.Label();
            this.label34 = new System.Windows.Forms.Label();
            this.label35 = new System.Windows.Forms.Label();
            this.label36 = new System.Windows.Forms.Label();
            this.labelIMUTempVal4 = new System.Windows.Forms.Label();
            this.labelAccZVal4 = new System.Windows.Forms.Label();
            this.labelAccYVal4 = new System.Windows.Forms.Label();
            this.labelAccXVal4 = new System.Windows.Forms.Label();
            this.labelMagZVal4 = new System.Windows.Forms.Label();
            this.labelMagYVal4 = new System.Windows.Forms.Label();
            this.labelMagXVal4 = new System.Windows.Forms.Label();
            this.labelGyroZVal4 = new System.Windows.Forms.Label();
            this.labelGyroYVal4 = new System.Windows.Forms.Label();
            this.labelGyroXVal4 = new System.Windows.Forms.Label();
            this.successByteLabel = new System.Windows.Forms.Label();
            this.label9 = new System.Windows.Forms.Label();
            this.groupBoxWarnings.SuspendLayout();
            this.groupBoxCRC.SuspendLayout();
            this.panel1.SuspendLayout();
            this.panel2.SuspendLayout();
            this.panel3.SuspendLayout();
            this.panel4.SuspendLayout();
            this.SuspendLayout();
            // 
            // serialComboBox
            // 
            this.serialComboBox.FormattingEnabled = true;
            this.serialComboBox.Items.AddRange(new object[] {
            "COM1",
            "COM2",
            "COM3",
            "COM4",
            "COM5",
            "COM6",
            "COM7",
            "COM8",
            "COM9",
            "COM10",
            "COM11",
            "COM12",
            "COM13",
            "COM14",
            "COM15",
            "COM16"});
            this.serialComboBox.Location = new System.Drawing.Point(10, 225);
            this.serialComboBox.MaxDropDownItems = 16;
            this.serialComboBox.Name = "serialComboBox";
            this.serialComboBox.Size = new System.Drawing.Size(129, 21);
            this.serialComboBox.TabIndex = 0;
            // 
            // labelSerialPort
            // 
            this.labelSerialPort.AutoSize = true;
            this.labelSerialPort.Font = new System.Drawing.Font("Arial", 12F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.labelSerialPort.Location = new System.Drawing.Point(28, 200);
            this.labelSerialPort.Name = "labelSerialPort";
            this.labelSerialPort.Size = new System.Drawing.Size(88, 19);
            this.labelSerialPort.TabIndex = 1;
            this.labelSerialPort.Text = "Serial Port";
            // 
            // openButton
            // 
            this.openButton.Font = new System.Drawing.Font("Arial", 14F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.openButton.Location = new System.Drawing.Point(10, 253);
            this.openButton.Name = "openButton";
            this.openButton.Size = new System.Drawing.Size(129, 39);
            this.openButton.TabIndex = 2;
            this.openButton.Text = "Open Port";
            this.openButton.UseVisualStyleBackColor = true;
            this.openButton.Click += new System.EventHandler(this.openButton_Click);
            // 
            // serialPort1
            // 
            this.serialPort1.BaudRate = 115200;
            this.serialPort1.PortName = "COM3";
            this.serialPort1.ReadBufferSize = 8192;
            // 
            // textBoxDebug
            // 
            this.textBoxDebug.Location = new System.Drawing.Point(10, 334);
            this.textBoxDebug.Multiline = true;
            this.textBoxDebug.Name = "textBoxDebug";
            this.textBoxDebug.ReadOnly = true;
            this.textBoxDebug.ScrollBars = System.Windows.Forms.ScrollBars.Vertical;
            this.textBoxDebug.Size = new System.Drawing.Size(275, 104);
            this.textBoxDebug.TabIndex = 7;
            // 
            // labelDebug
            // 
            this.labelDebug.AutoSize = true;
            this.labelDebug.Font = new System.Drawing.Font("Arial", 12F, System.Drawing.FontStyle.Bold, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
            this.labelDebug.Location = new System.Drawing.Point(10, 312);
            this.labelDebug.Name = "labelDebug";
            this.labelDebug.Size = new System.Drawing.Size(60, 19);
            this.labelDebug.TabIndex = 13;
            this.labelDebug.Text = "Debug";
            // 
            // buttonClearDebug
            // 
            this.buttonClearDebug.Location = new System.Drawing.Point(291, 334);
            this.buttonClearDebug.Name = "buttonClearDebug";
            this.buttonClearDebug.Size = new System.Drawing.Size(42, 22);
            this.buttonClearDebug.TabIndex = 21;
            this.buttonClearDebug.Text = "Clear";
            this.buttonClearDebug.UseVisualStyleBackColor = true;
            this.buttonClearDebug.Click += new System.EventHandler(this.buttonClearDebug_Click);
            // 
            // textBoxWarnings
            // 
            this.textBoxWarnings.Location = new System.Drawing.Point(0, 13);
            this.textBoxWarnings.Multiline = true;
            this.textBoxWarnings.Name = "textBoxWarnings";
            this.textBoxWarnings.ReadOnly = true;
            this.textBoxWarnings.ScrollBars = System.Windows.Forms.ScrollBars.Vertical;
            this.textBoxWarnings.Size = new System.Drawing.Size(262, 50);
            this.textBoxWarnings.TabIndex = 32;
            // 
            // groupBoxWarnings
            // 
            this.groupBoxWarnings.Controls.Add(this.textBoxWarnings);
            this.groupBoxWarnings.Location = new System.Drawing.Point(307, 228);
            this.groupBoxWarnings.Name = "groupBoxWarnings";
            this.groupBoxWarnings.Size = new System.Drawing.Size(262, 66);
            this.groupBoxWarnings.TabIndex = 33;
            this.groupBoxWarnings.TabStop = false;
            this.groupBoxWarnings.Text = "Warnings/Errors";
            // 
            // groupBoxCRC
            // 
            this.groupBoxCRC.Controls.Add(this.buttonResetCounters);
            this.groupBoxCRC.Controls.Add(this.labelEscCharFailsVal);
            this.groupBoxCRC.Controls.Add(this.labelCRCFailsVal);
            this.groupBoxCRC.Controls.Add(this.labelEscCharFails);
            this.groupBoxCRC.Controls.Add(this.labelCRCFails);
            this.groupBoxCRC.Location = new System.Drawing.Point(160, 200);
            this.groupBoxCRC.Name = "groupBoxCRC";
            this.groupBoxCRC.Size = new System.Drawing.Size(121, 117);
            this.groupBoxCRC.TabIndex = 36;
            this.groupBoxCRC.TabStop = false;
            this.groupBoxCRC.Text = "CRC Information";
            // 
            // buttonResetCounters
            // 
            this.buttonResetCounters.Location = new System.Drawing.Point(34, 74);
            this.buttonResetCounters.Name = "buttonResetCounters";
            this.buttonResetCounters.Size = new System.Drawing.Size(55, 35);
            this.buttonResetCounters.TabIndex = 39;
            this.buttonResetCounters.Text = "Reset";
            this.buttonResetCounters.UseVisualStyleBackColor = true;
            this.buttonResetCounters.Click += new System.EventHandler(this.buttonResetCounters_Click);
            // 
            // labelEscCharFailsVal
            // 
            this.labelEscCharFailsVal.AutoSize = true;
            this.labelEscCharFailsVal.Location = new System.Drawing.Point(76, 48);
            this.labelEscCharFailsVal.Name = "labelEscCharFailsVal";
            this.labelEscCharFailsVal.Size = new System.Drawing.Size(13, 13);
            this.labelEscCharFailsVal.TabIndex = 108;
            this.labelEscCharFailsVal.Text = "0";
            // 
            // labelCRCFailsVal
            // 
            this.labelCRCFailsVal.AutoSize = true;
            this.labelCRCFailsVal.Location = new System.Drawing.Point(76, 28);
            this.labelCRCFailsVal.Name = "labelCRCFailsVal";
            this.labelCRCFailsVal.Size = new System.Drawing.Size(13, 13);
            this.labelCRCFailsVal.TabIndex = 37;
            this.labelCRCFailsVal.Text = "0";
            // 
            // labelEscCharFails
            // 
            this.labelEscCharFails.AutoSize = true;
            this.labelEscCharFails.Location = new System.Drawing.Point(6, 48);
            this.labelEscCharFails.Name = "labelEscCharFails";
            this.labelEscCharFails.Size = new System.Drawing.Size(52, 13);
            this.labelEscCharFails.TabIndex = 107;
            this.labelEscCharFails.Text = "ESC Fails";
            // 
            // labelCRCFails
            // 
            this.labelCRCFails.AutoSize = true;
            this.labelCRCFails.Location = new System.Drawing.Point(6, 28);
            this.labelCRCFails.Name = "labelCRCFails";
            this.labelCRCFails.Size = new System.Drawing.Size(53, 13);
            this.labelCRCFails.TabIndex = 106;
            this.labelCRCFails.Text = "CRC Fails";
            // 
            // checkBoxLogging
            // 
            this.checkBoxLogging.AutoSize = true;
            this.checkBoxLogging.Location = new System.Drawing.Point(307, 203);
            this.checkBoxLogging.Name = "checkBoxLogging";
            this.checkBoxLogging.Size = new System.Drawing.Size(40, 17);
            this.checkBoxLogging.TabIndex = 37;
            this.checkBoxLogging.Text = "log";
            this.checkBoxLogging.UseVisualStyleBackColor = true;
            this.checkBoxLogging.CheckedChanged += new System.EventHandler(this.checkBoxLogging_CheckedChanged);
            // 
            // labelGyroXVal
            // 
            this.labelGyroXVal.AutoSize = true;
            this.labelGyroXVal.Location = new System.Drawing.Point(88, 17);
            this.labelGyroXVal.Name = "labelGyroXVal";
            this.labelGyroXVal.Size = new System.Drawing.Size(13, 13);
            this.labelGyroXVal.TabIndex = 38;
            this.labelGyroXVal.Text = "0";
            // 
            // labelGyroYVal
            // 
            this.labelGyroYVal.AutoSize = true;
            this.labelGyroYVal.Location = new System.Drawing.Point(88, 30);
            this.labelGyroYVal.Name = "labelGyroYVal";
            this.labelGyroYVal.Size = new System.Drawing.Size(13, 13);
            this.labelGyroYVal.TabIndex = 39;
            this.labelGyroYVal.Text = "0";
            // 
            // labelGyroZVal
            // 
            this.labelGyroZVal.AutoSize = true;
            this.labelGyroZVal.Location = new System.Drawing.Point(88, 43);
            this.labelGyroZVal.Name = "labelGyroZVal";
            this.labelGyroZVal.Size = new System.Drawing.Size(13, 13);
            this.labelGyroZVal.TabIndex = 40;
            this.labelGyroZVal.Text = "0";
            // 
            // labelMagXVal
            // 
            this.labelMagXVal.AutoSize = true;
            this.labelMagXVal.Location = new System.Drawing.Point(88, 56);
            this.labelMagXVal.Name = "labelMagXVal";
            this.labelMagXVal.Size = new System.Drawing.Size(13, 13);
            this.labelMagXVal.TabIndex = 41;
            this.labelMagXVal.Text = "0";
            // 
            // labelMagYVal
            // 
            this.labelMagYVal.AutoSize = true;
            this.labelMagYVal.Location = new System.Drawing.Point(88, 69);
            this.labelMagYVal.Name = "labelMagYVal";
            this.labelMagYVal.Size = new System.Drawing.Size(13, 13);
            this.labelMagYVal.TabIndex = 42;
            this.labelMagYVal.Text = "0";
            // 
            // labelMagZVal
            // 
            this.labelMagZVal.AutoSize = true;
            this.labelMagZVal.Location = new System.Drawing.Point(88, 82);
            this.labelMagZVal.Name = "labelMagZVal";
            this.labelMagZVal.Size = new System.Drawing.Size(13, 13);
            this.labelMagZVal.TabIndex = 43;
            this.labelMagZVal.Text = "0";
            // 
            // labelAccXVal
            // 
            this.labelAccXVal.AutoSize = true;
            this.labelAccXVal.Location = new System.Drawing.Point(88, 95);
            this.labelAccXVal.Name = "labelAccXVal";
            this.labelAccXVal.Size = new System.Drawing.Size(13, 13);
            this.labelAccXVal.TabIndex = 44;
            this.labelAccXVal.Text = "0";
            // 
            // labelAccYVal
            // 
            this.labelAccYVal.AutoSize = true;
            this.labelAccYVal.Location = new System.Drawing.Point(88, 108);
            this.labelAccYVal.Name = "labelAccYVal";
            this.labelAccYVal.Size = new System.Drawing.Size(13, 13);
            this.labelAccYVal.TabIndex = 45;
            this.labelAccYVal.Text = "0";
            // 
            // labelAccZVal
            // 
            this.labelAccZVal.AutoSize = true;
            this.labelAccZVal.Location = new System.Drawing.Point(88, 121);
            this.labelAccZVal.Name = "labelAccZVal";
            this.labelAccZVal.Size = new System.Drawing.Size(13, 13);
            this.labelAccZVal.TabIndex = 46;
            this.labelAccZVal.Text = "0";
            // 
            // labelIMUTempVal
            // 
            this.labelIMUTempVal.AutoSize = true;
            this.labelIMUTempVal.Location = new System.Drawing.Point(88, 134);
            this.labelIMUTempVal.Name = "labelIMUTempVal";
            this.labelIMUTempVal.Size = new System.Drawing.Size(13, 13);
            this.labelIMUTempVal.TabIndex = 47;
            this.labelIMUTempVal.Text = "0";
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(28, 17);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(29, 13);
            this.label1.TabIndex = 48;
            this.label1.Text = "Gyro";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(28, 56);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(28, 13);
            this.label2.TabIndex = 49;
            this.label2.Text = "Mag";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(28, 95);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(26, 13);
            this.label3.TabIndex = 50;
            this.label3.Text = "Acc";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(28, 134);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(34, 13);
            this.label4.TabIndex = 51;
            this.label4.Text = "Temp";
            // 
            // panel1
            // 
            this.panel1.BackColor = System.Drawing.SystemColors.ButtonHighlight;
            this.panel1.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.panel1.Controls.Add(this.ADCtempLABEL1);
            this.panel1.Controls.Add(this.label4);
            this.panel1.Controls.Add(this.label3);
            this.panel1.Controls.Add(this.label2);
            this.panel1.Controls.Add(this.label1);
            this.panel1.Controls.Add(this.labelIMUTempVal);
            this.panel1.Controls.Add(this.labelAccZVal);
            this.panel1.Controls.Add(this.labelAccYVal);
            this.panel1.Controls.Add(this.labelAccXVal);
            this.panel1.Controls.Add(this.labelMagZVal);
            this.panel1.Controls.Add(this.labelMagYVal);
            this.panel1.Controls.Add(this.labelMagXVal);
            this.panel1.Controls.Add(this.labelGyroZVal);
            this.panel1.Controls.Add(this.labelGyroYVal);
            this.panel1.Controls.Add(this.labelGyroXVal);
            this.panel1.Location = new System.Drawing.Point(14, 12);
            this.panel1.Name = "panel1";
            this.panel1.Size = new System.Drawing.Size(121, 174);
            this.panel1.TabIndex = 52;
            // 
            // ADCtempLABEL1
            // 
            this.ADCtempLABEL1.AutoSize = true;
            this.ADCtempLABEL1.Location = new System.Drawing.Point(88, 147);
            this.ADCtempLABEL1.Name = "ADCtempLABEL1";
            this.ADCtempLABEL1.Size = new System.Drawing.Size(13, 13);
            this.ADCtempLABEL1.TabIndex = 52;
            this.ADCtempLABEL1.Text = "0";
            // 
            // panel2
            // 
            this.panel2.BackColor = System.Drawing.SystemColors.ButtonHighlight;
            this.panel2.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.panel2.Controls.Add(this.ADCtempLABEL2);
            this.panel2.Controls.Add(this.label5);
            this.panel2.Controls.Add(this.label6);
            this.panel2.Controls.Add(this.label7);
            this.panel2.Controls.Add(this.label8);
            this.panel2.Controls.Add(this.labelIMUTempVal2);
            this.panel2.Controls.Add(this.labelAccZVal2);
            this.panel2.Controls.Add(this.labelAccYVal2);
            this.panel2.Controls.Add(this.labelAccXVal2);
            this.panel2.Controls.Add(this.labelMagZVal2);
            this.panel2.Controls.Add(this.labelMagYVal2);
            this.panel2.Controls.Add(this.labelMagXVal2);
            this.panel2.Controls.Add(this.labelGyroZVal2);
            this.panel2.Controls.Add(this.labelGyroYVal2);
            this.panel2.Controls.Add(this.labelGyroXVal2);
            this.panel2.Location = new System.Drawing.Point(160, 12);
            this.panel2.Name = "panel2";
            this.panel2.Size = new System.Drawing.Size(121, 174);
            this.panel2.TabIndex = 53;
            // 
            // ADCtempLABEL2
            // 
            this.ADCtempLABEL2.AutoSize = true;
            this.ADCtempLABEL2.Location = new System.Drawing.Point(88, 147);
            this.ADCtempLABEL2.Name = "ADCtempLABEL2";
            this.ADCtempLABEL2.Size = new System.Drawing.Size(13, 13);
            this.ADCtempLABEL2.TabIndex = 52;
            this.ADCtempLABEL2.Text = "0";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(28, 134);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(34, 13);
            this.label5.TabIndex = 51;
            this.label5.Text = "Temp";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(28, 95);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(26, 13);
            this.label6.TabIndex = 50;
            this.label6.Text = "Acc";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(28, 56);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(28, 13);
            this.label7.TabIndex = 49;
            this.label7.Text = "Mag";
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(28, 17);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(29, 13);
            this.label8.TabIndex = 48;
            this.label8.Text = "Gyro";
            // 
            // labelIMUTempVal2
            // 
            this.labelIMUTempVal2.AutoSize = true;
            this.labelIMUTempVal2.Location = new System.Drawing.Point(88, 134);
            this.labelIMUTempVal2.Name = "labelIMUTempVal2";
            this.labelIMUTempVal2.Size = new System.Drawing.Size(13, 13);
            this.labelIMUTempVal2.TabIndex = 47;
            this.labelIMUTempVal2.Text = "0";
            // 
            // labelAccZVal2
            // 
            this.labelAccZVal2.AutoSize = true;
            this.labelAccZVal2.Location = new System.Drawing.Point(88, 121);
            this.labelAccZVal2.Name = "labelAccZVal2";
            this.labelAccZVal2.Size = new System.Drawing.Size(13, 13);
            this.labelAccZVal2.TabIndex = 46;
            this.labelAccZVal2.Text = "0";
            // 
            // labelAccYVal2
            // 
            this.labelAccYVal2.AutoSize = true;
            this.labelAccYVal2.Location = new System.Drawing.Point(88, 108);
            this.labelAccYVal2.Name = "labelAccYVal2";
            this.labelAccYVal2.Size = new System.Drawing.Size(13, 13);
            this.labelAccYVal2.TabIndex = 45;
            this.labelAccYVal2.Text = "0";
            // 
            // labelAccXVal2
            // 
            this.labelAccXVal2.AutoSize = true;
            this.labelAccXVal2.Location = new System.Drawing.Point(88, 95);
            this.labelAccXVal2.Name = "labelAccXVal2";
            this.labelAccXVal2.Size = new System.Drawing.Size(13, 13);
            this.labelAccXVal2.TabIndex = 44;
            this.labelAccXVal2.Text = "0";
            // 
            // labelMagZVal2
            // 
            this.labelMagZVal2.AutoSize = true;
            this.labelMagZVal2.Location = new System.Drawing.Point(88, 82);
            this.labelMagZVal2.Name = "labelMagZVal2";
            this.labelMagZVal2.Size = new System.Drawing.Size(13, 13);
            this.labelMagZVal2.TabIndex = 43;
            this.labelMagZVal2.Text = "0";
            // 
            // labelMagYVal2
            // 
            this.labelMagYVal2.AutoSize = true;
            this.labelMagYVal2.Location = new System.Drawing.Point(88, 69);
            this.labelMagYVal2.Name = "labelMagYVal2";
            this.labelMagYVal2.Size = new System.Drawing.Size(13, 13);
            this.labelMagYVal2.TabIndex = 42;
            this.labelMagYVal2.Text = "0";
            // 
            // labelMagXVal2
            // 
            this.labelMagXVal2.AutoSize = true;
            this.labelMagXVal2.Location = new System.Drawing.Point(88, 56);
            this.labelMagXVal2.Name = "labelMagXVal2";
            this.labelMagXVal2.Size = new System.Drawing.Size(13, 13);
            this.labelMagXVal2.TabIndex = 41;
            this.labelMagXVal2.Text = "0";
            // 
            // labelGyroZVal2
            // 
            this.labelGyroZVal2.AutoSize = true;
            this.labelGyroZVal2.Location = new System.Drawing.Point(88, 43);
            this.labelGyroZVal2.Name = "labelGyroZVal2";
            this.labelGyroZVal2.Size = new System.Drawing.Size(13, 13);
            this.labelGyroZVal2.TabIndex = 40;
            this.labelGyroZVal2.Text = "0";
            // 
            // labelGyroYVal2
            // 
            this.labelGyroYVal2.AutoSize = true;
            this.labelGyroYVal2.Location = new System.Drawing.Point(88, 30);
            this.labelGyroYVal2.Name = "labelGyroYVal2";
            this.labelGyroYVal2.Size = new System.Drawing.Size(13, 13);
            this.labelGyroYVal2.TabIndex = 39;
            this.labelGyroYVal2.Text = "0";
            // 
            // labelGyroXVal2
            // 
            this.labelGyroXVal2.AutoSize = true;
            this.labelGyroXVal2.Location = new System.Drawing.Point(88, 17);
            this.labelGyroXVal2.Name = "labelGyroXVal2";
            this.labelGyroXVal2.Size = new System.Drawing.Size(13, 13);
            this.labelGyroXVal2.TabIndex = 38;
            this.labelGyroXVal2.Text = "0";
            // 
            // panel3
            // 
            this.panel3.BackColor = System.Drawing.SystemColors.ButtonHighlight;
            this.panel3.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.panel3.Controls.Add(this.ADCtempLABEL3);
            this.panel3.Controls.Add(this.label19);
            this.panel3.Controls.Add(this.label20);
            this.panel3.Controls.Add(this.label21);
            this.panel3.Controls.Add(this.label22);
            this.panel3.Controls.Add(this.labelIMUTempVal3);
            this.panel3.Controls.Add(this.labelAccZVal3);
            this.panel3.Controls.Add(this.labelAccYVal3);
            this.panel3.Controls.Add(this.labelAccXVal3);
            this.panel3.Controls.Add(this.labelMagZVal3);
            this.panel3.Controls.Add(this.labelMagYVal3);
            this.panel3.Controls.Add(this.labelMagXVal3);
            this.panel3.Controls.Add(this.labelGyroZVal3);
            this.panel3.Controls.Add(this.labelGyroYVal3);
            this.panel3.Controls.Add(this.labelGyroXVal3);
            this.panel3.Location = new System.Drawing.Point(307, 12);
            this.panel3.Name = "panel3";
            this.panel3.Size = new System.Drawing.Size(121, 174);
            this.panel3.TabIndex = 53;
            // 
            // ADCtempLABEL3
            // 
            this.ADCtempLABEL3.AutoSize = true;
            this.ADCtempLABEL3.Location = new System.Drawing.Point(88, 147);
            this.ADCtempLABEL3.Name = "ADCtempLABEL3";
            this.ADCtempLABEL3.Size = new System.Drawing.Size(13, 13);
            this.ADCtempLABEL3.TabIndex = 52;
            this.ADCtempLABEL3.Text = "0";
            // 
            // label19
            // 
            this.label19.AutoSize = true;
            this.label19.Location = new System.Drawing.Point(28, 134);
            this.label19.Name = "label19";
            this.label19.Size = new System.Drawing.Size(34, 13);
            this.label19.TabIndex = 51;
            this.label19.Text = "Temp";
            // 
            // label20
            // 
            this.label20.AutoSize = true;
            this.label20.Location = new System.Drawing.Point(28, 95);
            this.label20.Name = "label20";
            this.label20.Size = new System.Drawing.Size(26, 13);
            this.label20.TabIndex = 50;
            this.label20.Text = "Acc";
            // 
            // label21
            // 
            this.label21.AutoSize = true;
            this.label21.Location = new System.Drawing.Point(28, 56);
            this.label21.Name = "label21";
            this.label21.Size = new System.Drawing.Size(28, 13);
            this.label21.TabIndex = 49;
            this.label21.Text = "Mag";
            // 
            // label22
            // 
            this.label22.AutoSize = true;
            this.label22.Location = new System.Drawing.Point(28, 17);
            this.label22.Name = "label22";
            this.label22.Size = new System.Drawing.Size(29, 13);
            this.label22.TabIndex = 48;
            this.label22.Text = "Gyro";
            // 
            // labelIMUTempVal3
            // 
            this.labelIMUTempVal3.AutoSize = true;
            this.labelIMUTempVal3.Location = new System.Drawing.Point(88, 134);
            this.labelIMUTempVal3.Name = "labelIMUTempVal3";
            this.labelIMUTempVal3.Size = new System.Drawing.Size(13, 13);
            this.labelIMUTempVal3.TabIndex = 47;
            this.labelIMUTempVal3.Text = "0";
            // 
            // labelAccZVal3
            // 
            this.labelAccZVal3.AutoSize = true;
            this.labelAccZVal3.Location = new System.Drawing.Point(88, 121);
            this.labelAccZVal3.Name = "labelAccZVal3";
            this.labelAccZVal3.Size = new System.Drawing.Size(13, 13);
            this.labelAccZVal3.TabIndex = 46;
            this.labelAccZVal3.Text = "0";
            // 
            // labelAccYVal3
            // 
            this.labelAccYVal3.AutoSize = true;
            this.labelAccYVal3.Location = new System.Drawing.Point(88, 108);
            this.labelAccYVal3.Name = "labelAccYVal3";
            this.labelAccYVal3.Size = new System.Drawing.Size(13, 13);
            this.labelAccYVal3.TabIndex = 45;
            this.labelAccYVal3.Text = "0";
            // 
            // labelAccXVal3
            // 
            this.labelAccXVal3.AutoSize = true;
            this.labelAccXVal3.Location = new System.Drawing.Point(88, 95);
            this.labelAccXVal3.Name = "labelAccXVal3";
            this.labelAccXVal3.Size = new System.Drawing.Size(13, 13);
            this.labelAccXVal3.TabIndex = 44;
            this.labelAccXVal3.Text = "0";
            // 
            // labelMagZVal3
            // 
            this.labelMagZVal3.AutoSize = true;
            this.labelMagZVal3.Location = new System.Drawing.Point(88, 82);
            this.labelMagZVal3.Name = "labelMagZVal3";
            this.labelMagZVal3.Size = new System.Drawing.Size(13, 13);
            this.labelMagZVal3.TabIndex = 43;
            this.labelMagZVal3.Text = "0";
            // 
            // labelMagYVal3
            // 
            this.labelMagYVal3.AutoSize = true;
            this.labelMagYVal3.Location = new System.Drawing.Point(88, 69);
            this.labelMagYVal3.Name = "labelMagYVal3";
            this.labelMagYVal3.Size = new System.Drawing.Size(13, 13);
            this.labelMagYVal3.TabIndex = 42;
            this.labelMagYVal3.Text = "0";
            // 
            // labelMagXVal3
            // 
            this.labelMagXVal3.AutoSize = true;
            this.labelMagXVal3.Location = new System.Drawing.Point(88, 56);
            this.labelMagXVal3.Name = "labelMagXVal3";
            this.labelMagXVal3.Size = new System.Drawing.Size(13, 13);
            this.labelMagXVal3.TabIndex = 41;
            this.labelMagXVal3.Text = "0";
            // 
            // labelGyroZVal3
            // 
            this.labelGyroZVal3.AutoSize = true;
            this.labelGyroZVal3.Location = new System.Drawing.Point(88, 43);
            this.labelGyroZVal3.Name = "labelGyroZVal3";
            this.labelGyroZVal3.Size = new System.Drawing.Size(13, 13);
            this.labelGyroZVal3.TabIndex = 40;
            this.labelGyroZVal3.Text = "0";
            // 
            // labelGyroYVal3
            // 
            this.labelGyroYVal3.AutoSize = true;
            this.labelGyroYVal3.Location = new System.Drawing.Point(88, 30);
            this.labelGyroYVal3.Name = "labelGyroYVal3";
            this.labelGyroYVal3.Size = new System.Drawing.Size(13, 13);
            this.labelGyroYVal3.TabIndex = 39;
            this.labelGyroYVal3.Text = "0";
            // 
            // labelGyroXVal3
            // 
            this.labelGyroXVal3.AutoSize = true;
            this.labelGyroXVal3.Location = new System.Drawing.Point(88, 17);
            this.labelGyroXVal3.Name = "labelGyroXVal3";
            this.labelGyroXVal3.Size = new System.Drawing.Size(13, 13);
            this.labelGyroXVal3.TabIndex = 38;
            this.labelGyroXVal3.Text = "0";
            // 
            // panel4
            // 
            this.panel4.BackColor = System.Drawing.SystemColors.ButtonHighlight;
            this.panel4.BorderStyle = System.Windows.Forms.BorderStyle.FixedSingle;
            this.panel4.Controls.Add(this.ADCtempLABEL4);
            this.panel4.Controls.Add(this.label33);
            this.panel4.Controls.Add(this.label34);
            this.panel4.Controls.Add(this.label35);
            this.panel4.Controls.Add(this.label36);
            this.panel4.Controls.Add(this.labelIMUTempVal4);
            this.panel4.Controls.Add(this.labelAccZVal4);
            this.panel4.Controls.Add(this.labelAccYVal4);
            this.panel4.Controls.Add(this.labelAccXVal4);
            this.panel4.Controls.Add(this.labelMagZVal4);
            this.panel4.Controls.Add(this.labelMagYVal4);
            this.panel4.Controls.Add(this.labelMagXVal4);
            this.panel4.Controls.Add(this.labelGyroZVal4);
            this.panel4.Controls.Add(this.labelGyroYVal4);
            this.panel4.Controls.Add(this.labelGyroXVal4);
            this.panel4.Location = new System.Drawing.Point(448, 12);
            this.panel4.Name = "panel4";
            this.panel4.Size = new System.Drawing.Size(121, 174);
            this.panel4.TabIndex = 53;
            // 
            // ADCtempLABEL4
            // 
            this.ADCtempLABEL4.AutoSize = true;
            this.ADCtempLABEL4.Location = new System.Drawing.Point(88, 147);
            this.ADCtempLABEL4.Name = "ADCtempLABEL4";
            this.ADCtempLABEL4.Size = new System.Drawing.Size(13, 13);
            this.ADCtempLABEL4.TabIndex = 53;
            this.ADCtempLABEL4.Text = "0";
            // 
            // label33
            // 
            this.label33.AutoSize = true;
            this.label33.Location = new System.Drawing.Point(28, 134);
            this.label33.Name = "label33";
            this.label33.Size = new System.Drawing.Size(34, 13);
            this.label33.TabIndex = 51;
            this.label33.Text = "Temp";
            // 
            // label34
            // 
            this.label34.AutoSize = true;
            this.label34.Location = new System.Drawing.Point(28, 95);
            this.label34.Name = "label34";
            this.label34.Size = new System.Drawing.Size(26, 13);
            this.label34.TabIndex = 50;
            this.label34.Text = "Acc";
            // 
            // label35
            // 
            this.label35.AutoSize = true;
            this.label35.Location = new System.Drawing.Point(28, 56);
            this.label35.Name = "label35";
            this.label35.Size = new System.Drawing.Size(28, 13);
            this.label35.TabIndex = 49;
            this.label35.Text = "Mag";
            // 
            // label36
            // 
            this.label36.AutoSize = true;
            this.label36.Location = new System.Drawing.Point(28, 17);
            this.label36.Name = "label36";
            this.label36.Size = new System.Drawing.Size(29, 13);
            this.label36.TabIndex = 48;
            this.label36.Text = "Gyro";
            // 
            // labelIMUTempVal4
            // 
            this.labelIMUTempVal4.AutoSize = true;
            this.labelIMUTempVal4.Location = new System.Drawing.Point(88, 134);
            this.labelIMUTempVal4.Name = "labelIMUTempVal4";
            this.labelIMUTempVal4.Size = new System.Drawing.Size(13, 13);
            this.labelIMUTempVal4.TabIndex = 47;
            this.labelIMUTempVal4.Text = "0";
            // 
            // labelAccZVal4
            // 
            this.labelAccZVal4.AutoSize = true;
            this.labelAccZVal4.Location = new System.Drawing.Point(88, 121);
            this.labelAccZVal4.Name = "labelAccZVal4";
            this.labelAccZVal4.Size = new System.Drawing.Size(13, 13);
            this.labelAccZVal4.TabIndex = 46;
            this.labelAccZVal4.Text = "0";
            // 
            // labelAccYVal4
            // 
            this.labelAccYVal4.AutoSize = true;
            this.labelAccYVal4.Location = new System.Drawing.Point(88, 108);
            this.labelAccYVal4.Name = "labelAccYVal4";
            this.labelAccYVal4.Size = new System.Drawing.Size(13, 13);
            this.labelAccYVal4.TabIndex = 45;
            this.labelAccYVal4.Text = "0";
            // 
            // labelAccXVal4
            // 
            this.labelAccXVal4.AutoSize = true;
            this.labelAccXVal4.Location = new System.Drawing.Point(88, 95);
            this.labelAccXVal4.Name = "labelAccXVal4";
            this.labelAccXVal4.Size = new System.Drawing.Size(13, 13);
            this.labelAccXVal4.TabIndex = 44;
            this.labelAccXVal4.Text = "0";
            // 
            // labelMagZVal4
            // 
            this.labelMagZVal4.AutoSize = true;
            this.labelMagZVal4.Location = new System.Drawing.Point(88, 82);
            this.labelMagZVal4.Name = "labelMagZVal4";
            this.labelMagZVal4.Size = new System.Drawing.Size(13, 13);
            this.labelMagZVal4.TabIndex = 43;
            this.labelMagZVal4.Text = "0";
            // 
            // labelMagYVal4
            // 
            this.labelMagYVal4.AutoSize = true;
            this.labelMagYVal4.Location = new System.Drawing.Point(88, 69);
            this.labelMagYVal4.Name = "labelMagYVal4";
            this.labelMagYVal4.Size = new System.Drawing.Size(13, 13);
            this.labelMagYVal4.TabIndex = 42;
            this.labelMagYVal4.Text = "0";
            // 
            // labelMagXVal4
            // 
            this.labelMagXVal4.AutoSize = true;
            this.labelMagXVal4.Location = new System.Drawing.Point(88, 56);
            this.labelMagXVal4.Name = "labelMagXVal4";
            this.labelMagXVal4.Size = new System.Drawing.Size(13, 13);
            this.labelMagXVal4.TabIndex = 41;
            this.labelMagXVal4.Text = "0";
            // 
            // labelGyroZVal4
            // 
            this.labelGyroZVal4.AutoSize = true;
            this.labelGyroZVal4.Location = new System.Drawing.Point(88, 43);
            this.labelGyroZVal4.Name = "labelGyroZVal4";
            this.labelGyroZVal4.Size = new System.Drawing.Size(13, 13);
            this.labelGyroZVal4.TabIndex = 40;
            this.labelGyroZVal4.Text = "0";
            // 
            // labelGyroYVal4
            // 
            this.labelGyroYVal4.AutoSize = true;
            this.labelGyroYVal4.Location = new System.Drawing.Point(88, 30);
            this.labelGyroYVal4.Name = "labelGyroYVal4";
            this.labelGyroYVal4.Size = new System.Drawing.Size(13, 13);
            this.labelGyroYVal4.TabIndex = 39;
            this.labelGyroYVal4.Text = "0";
            // 
            // labelGyroXVal4
            // 
            this.labelGyroXVal4.AutoSize = true;
            this.labelGyroXVal4.Location = new System.Drawing.Point(88, 17);
            this.labelGyroXVal4.Name = "labelGyroXVal4";
            this.labelGyroXVal4.Size = new System.Drawing.Size(13, 13);
            this.labelGyroXVal4.TabIndex = 38;
            this.labelGyroXVal4.Text = "0";
            // 
            // successByteLabel
            // 
            this.successByteLabel.AutoSize = true;
            this.successByteLabel.Location = new System.Drawing.Point(542, 200);
            this.successByteLabel.Name = "successByteLabel";
            this.successByteLabel.Size = new System.Drawing.Size(13, 13);
            this.successByteLabel.TabIndex = 54;
            this.successByteLabel.Text = "0";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(454, 200);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(82, 13);
            this.label9.TabIndex = 55;
            this.label9.Text = "packet success";
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(606, 445);
            this.Controls.Add(this.label9);
            this.Controls.Add(this.successByteLabel);
            this.Controls.Add(this.panel4);
            this.Controls.Add(this.panel3);
            this.Controls.Add(this.panel2);
            this.Controls.Add(this.panel1);
            this.Controls.Add(this.checkBoxLogging);
            this.Controls.Add(this.groupBoxCRC);
            this.Controls.Add(this.groupBoxWarnings);
            this.Controls.Add(this.buttonClearDebug);
            this.Controls.Add(this.labelDebug);
            this.Controls.Add(this.textBoxDebug);
            this.Controls.Add(this.openButton);
            this.Controls.Add(this.labelSerialPort);
            this.Controls.Add(this.serialComboBox);
            this.Name = "Form1";
            this.Text = "Cheetah Data Logger App";
            this.Load += new System.EventHandler(this.Form1_Load);
            this.groupBoxWarnings.ResumeLayout(false);
            this.groupBoxWarnings.PerformLayout();
            this.groupBoxCRC.ResumeLayout(false);
            this.groupBoxCRC.PerformLayout();
            this.panel1.ResumeLayout(false);
            this.panel1.PerformLayout();
            this.panel2.ResumeLayout(false);
            this.panel2.PerformLayout();
            this.panel3.ResumeLayout(false);
            this.panel3.PerformLayout();
            this.panel4.ResumeLayout(false);
            this.panel4.PerformLayout();
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Label labelSerialPort;
        private System.Windows.Forms.Button openButton;
        private System.IO.Ports.SerialPort serialPort1;
        private System.Windows.Forms.TextBox textBoxDebug;
        private System.Windows.Forms.Label labelDebug;
        private System.Windows.Forms.ComboBox serialComboBox;
        private System.Windows.Forms.Button buttonClearDebug;
        private System.Windows.Forms.TextBox textBoxWarnings;
        private System.Windows.Forms.GroupBox groupBoxWarnings;
        private System.Windows.Forms.GroupBox groupBoxCRC;
        private System.Windows.Forms.Label labelEscCharFailsVal;
        private System.Windows.Forms.Label labelCRCFailsVal;
        private System.Windows.Forms.Label labelEscCharFails;
        private System.Windows.Forms.Label labelCRCFails;
        private System.Windows.Forms.Button buttonResetCounters;
        private System.Windows.Forms.CheckBox checkBoxLogging;
        private System.Windows.Forms.Label labelGyroXVal;
        private System.Windows.Forms.Label labelGyroYVal;
        private System.Windows.Forms.Label labelGyroZVal;
        private System.Windows.Forms.Label labelMagXVal;
        private System.Windows.Forms.Label labelMagYVal;
        private System.Windows.Forms.Label labelMagZVal;
        private System.Windows.Forms.Label labelAccXVal;
        private System.Windows.Forms.Label labelAccYVal;
        private System.Windows.Forms.Label labelAccZVal;
        private System.Windows.Forms.Label labelIMUTempVal;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Panel panel1;
        private System.Windows.Forms.Panel panel2;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.Label labelIMUTempVal2;
        private System.Windows.Forms.Label labelAccZVal2;
        private System.Windows.Forms.Label labelAccYVal2;
        private System.Windows.Forms.Label labelAccXVal2;
        private System.Windows.Forms.Label labelMagZVal2;
        private System.Windows.Forms.Label labelMagYVal2;
        private System.Windows.Forms.Label labelMagXVal2;
        private System.Windows.Forms.Label labelGyroZVal2;
        private System.Windows.Forms.Label labelGyroYVal2;
        private System.Windows.Forms.Label labelGyroXVal2;
        private System.Windows.Forms.Panel panel3;
        private System.Windows.Forms.Label label19;
        private System.Windows.Forms.Label label20;
        private System.Windows.Forms.Label label21;
        private System.Windows.Forms.Label label22;
        private System.Windows.Forms.Label labelIMUTempVal3;
        private System.Windows.Forms.Label labelAccZVal3;
        private System.Windows.Forms.Label labelAccYVal3;
        private System.Windows.Forms.Label labelAccXVal3;
        private System.Windows.Forms.Label labelMagZVal3;
        private System.Windows.Forms.Label labelMagYVal3;
        private System.Windows.Forms.Label labelMagXVal3;
        private System.Windows.Forms.Label labelGyroZVal3;
        private System.Windows.Forms.Label labelGyroYVal3;
        private System.Windows.Forms.Label labelGyroXVal3;
        private System.Windows.Forms.Panel panel4;
        private System.Windows.Forms.Label label33;
        private System.Windows.Forms.Label label34;
        private System.Windows.Forms.Label label35;
        private System.Windows.Forms.Label label36;
        private System.Windows.Forms.Label labelIMUTempVal4;
        private System.Windows.Forms.Label labelAccZVal4;
        private System.Windows.Forms.Label labelAccYVal4;
        private System.Windows.Forms.Label labelAccXVal4;
        private System.Windows.Forms.Label labelMagZVal4;
        private System.Windows.Forms.Label labelMagYVal4;
        private System.Windows.Forms.Label labelMagXVal4;
        private System.Windows.Forms.Label labelGyroZVal4;
        private System.Windows.Forms.Label labelGyroYVal4;
        private System.Windows.Forms.Label labelGyroXVal4;
        private System.Windows.Forms.Label ADCtempLABEL1;
        private System.Windows.Forms.Label ADCtempLABEL2;
        private System.Windows.Forms.Label ADCtempLABEL3;
        private System.Windows.Forms.Label ADCtempLABEL4;
        private System.Windows.Forms.Label successByteLabel;
        private System.Windows.Forms.Label label9;
    }
}

