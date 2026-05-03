nRF24L01-Digital - Trasnsreciever 

A 2.4GHz digital voice transceiver using Arduino, nRF24L01+, and GFSK modulation, featuring MATLAB simulations for BT-product optimization and spectral efficiency.

Compenents Used:(For transmitter and reciever)
Arduino Nano / Uno	            2	      Core microcontroller for audio ADC sampling and SPI management.
nRF24L01+ Module    	          2	      2.4GHz RF Transceiver.
Electret Microphone Amplifier	  2	      (e.g., MAX4466) To capture voice and amplify it for the Arduino ADC.
PAM8403 Amplifier	              2      	Class-D audio amplifier to drive the speaker.
8Ω Speaker	                    2	      Audio output (0.5W to 3W).
Push Button	                    2	      Used for the Push-to-Talk (PTT) mechanism.
10µF - 100µF Capacitor	        2	      Bypass capacitor across nRF24L01+ power pins (CRITICAL for stability).
Resistors & Capacitors	        Set    	Used to build the custom RC Low-Pass Filter.

Repository Structure
├── README.md
├── Documentation/
│   ├── ckt_diagram.jpeg
│   ├── flowchart.png
│   ├── flowchart_2.jpeg

├── Hardware/
│   ├── Walkie_Talkie_Main.ino
│   ├── prototype_img_1.jpeg
│   ├── RF24Audio-1.0.zip
│   ├── Arduino_Walkie_Talkie_nRF24L01
|   |   ├── Arduino_Walkie_Talkie_nRF24L01.ino

├── matlab_sim/
│   ├── fsk_vs_gfsk
|   |   ├── Figure.png
|   |   ├── Figure_2.png
|   |   ├── fsk_vs_gfsk.m
│   ├── gfsk_sims
|   |   ├── ber.png
|   |   ├── psd.png
|   |   ├── Gfsk_sim.m
