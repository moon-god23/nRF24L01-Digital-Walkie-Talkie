#include <SPI.h>
#include <RF24.h>
#include <RF24Audio.h>
#include "printf.h"

RF24 radio(7, 8); // CE = 7, CSN = 8
RF24Audio rfAudio(radio, 0); 
const int talkButton = 3; // Button wired between Pin 3 and GND

// State tracking to prevent rapid switching
bool isTransmitting = false;

void setup() {      
  Serial.begin(115200);
  printf_begin();
  
  // Initialize hardware
  radio.begin();
  rfAudio.begin();

  // Network and Power Settings
  radio.setChannel(115); // Dodge Wi-Fi interference
  radio.setPALevel(RF24_PA_MAX); // Max power (requires capacitors)
  rfAudio.setVolume(4); // Cap volume to prevent PAM8403 speaker blowout
  
  // Initialize button with internal pull-up resistor
  pinMode(talkButton, INPUT_PULLUP);
  
  // Default state is receiving
  Serial.println("Walkie-Talkie Booted. Default: Listening...");
  rfAudio.receive();
}

void loop() {
  // LOW means the button is actively being pressed (grounded)
  bool buttonPressed = (digitalRead(talkButton) == LOW);

  // If button is held down AND we aren't already transmitting -> Switch to TX
  if (buttonPressed && !isTransmitting) {
    delay(30); // 30ms Debounce buffer to ignore physical switch bounce
    if (digitalRead(talkButton) == LOW) { 
      Serial.println("Transmitting...");
      rfAudio.transmit();
      isTransmitting = true;
    }
  }
  
  // If button is released AND we are currently transmitting -> Switch back to RX
  else if (!buttonPressed && isTransmitting) {
    delay(30); // 30ms Debounce buffer
    if (digitalRead(talkButton) == HIGH) { 
      Serial.println("Listening...");
      rfAudio.receive();
      isTransmitting = false;
    }
  }
}