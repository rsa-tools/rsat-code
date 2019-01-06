/*========================================================================

  File: chronometer.c

  Program REA - Version 1.1 - July 1999

  ========================================================================

  This module contains an implementation of a chronometer, designed to
  measure cumulated running time of seleted portions of code, with a
  precision of 0.01 seconds. 

  See the file chronometer.h for more information

  ========================================================================

    Copyright (C) 1999 Victor Jimenez and Andres Marzal

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program (file COPYING); if not, write to the Free
    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    You can contact the authors at:

           Victor Jimenez and Andres Marzal
           Dept. de Informatica, Universitat Jaume I, Castellon, Spain
           {vjimenez,amarzal}@inf.uji.es

    Updated information about this software will be available at the
    following www address:

           http://terra.act.uji.es/REA
 
  ========================================================================= */


#include <stdio.h>
#include <sys/times.h>  

#include <unistd.h>

#ifdef DEBUG
#include <assert.h>
#endif

#define CLOCK_STOPPED 0
#define CLOCK_RUNNING 1

/*============================================================================
  GLOBAL VARIABLES
  ============================================================================*/

int clockStatus = CLOCK_STOPPED;
struct tms clockTicks;
int totalClockTicks = 0;


/*============================================================================*/
void ClockReset () {
  /* Sets the cumulated time to 0 and starts cumulating time */

   totalClockTicks = 0;
   clockStatus = CLOCK_RUNNING;
   while (times(&clockTicks) == -1);
   
}


/*============================================================================*/
void ClockStop () {
  /* Stops cumulating time and keeps the currect value */

   struct tms clockTicks2;

   while (times(&clockTicks2) == -1);

#ifdef DEBUG   
   assert(clockStatus == CLOCK_RUNNING);
#endif

   totalClockTicks += (clockTicks2.tms_utime - clockTicks.tms_utime);
   clockTicks = clockTicks2;
   clockStatus = CLOCK_STOPPED;

}


/*============================================================================*/
void ClockContinue () {
  /* Continues cumulating time to the value when it was stopped */

#ifdef DEBUG
   assert (clockStatus == CLOCK_STOPPED);
#endif

   clockStatus = CLOCK_RUNNING;
   while (times(&clockTicks) == -1);

}


/*============================================================================*/
float ClockTotal () {
  /* Returns the cumulated time, in seconds */

   struct tms clockTicks2;

   while (times(&clockTicks2) == -1);

   if (clockStatus == CLOCK_RUNNING) {
      totalClockTicks += (clockTicks2.tms_utime - clockTicks.tms_utime);
      clockTicks = clockTicks2;
   }
   
   /* sysconf (_SC_CLK_TCK) returns the clicks per second */
   return (totalClockTicks / (float) sysconf (_SC_CLK_TCK));
   
}
