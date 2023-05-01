/*------ Display Connections ------
 * D = SDA = A4
 * C = SCL = A5
 * 
 */

#include <Wire.h> // Enable this line if using Arduino Uno, Mega, etc.
#include <Adafruit_GFX.h>
#include "Adafruit_LEDBackpack.h"

Adafruit_7segment matrix = Adafruit_7segment();

int TachPin = 4;
unsigned long T;    //stores time in [ms]
unsigned long oldT = 0;
int updateTime = 5000;  //time to wait while shit spins
float t = 0;
int revs = 0;
int switchOld = 1;
int slots = 1;        //number of slots on encoder wheel
float Ripumms;

void setup() {
  pinMode(TachPin, INPUT);
  //#ifndef __AVR_ATtiny85__
  //Serial.begin(9600);
  //Serial.println("7 Segment Backpack Test");
//#endif
  matrix.begin(0x70);
  matrix.print(10000);
  matrix.writeDisplay();
}

void loop() {
  T = millis();

  int SwitchState = digitalRead(TachPin);
  if(SwitchState != switchOld){
    revs++;   //This # is twice as high b/c there's 2 state changes
  }
  

  if(T-oldT>updateTime){
    t = (T-oldT)/1000;        //Time [s] since last update = updateTime
    // do shit here
    Ripumms = revs*60/t/slots/2;      // RPM's = (revolutions/ time [s])*(60 [s]/1 [min])*(1/# of encoder slots)*(1/2)
    /*Serial.print("Revs = ");          //    The 1/2 is from 2 state changes per slot
    Serial.print(revs);
    Serial.print("\t");
    Serial.print("Time = ");
    Serial.print(t,3);
    Serial.print("\t");
    Serial.print("RPM = ");
    Serial.println(Ripumms);
*/
    matrix.print(int(Ripumms));
    matrix.writeDisplay();
    
    revs = 0;
    oldT = T;
  }

  switchOld = SwitchState;
}
