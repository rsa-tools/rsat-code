/*========================================================================

  File: chronometer.h

  Program REA - Version 1.1 - July 1999

  ========================================================================

  This module contains an implementation of a chronometer, designed to
  measure cumulated running time of seleted portions of code, with a
  precision of 0.01 seconds. 

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
 
  =========================================================================

  This module includes the following functions:

     ClockReset() sets the cumulated time to 0 and starts cumulating time
     ClockStop() stops cumulating time and keeps the currect value
     ClockContinue() continues cumulating time to the value when it was stopped
     ClockTotal() returns the cumulated time, in seconds

  These functions use the standard C function times (that returns the number
  of clock ticks that have elapsed since the system has been up) and takes
  the CPU time dedicated to the process and its descendents, but not to the
  operating system.

  In Example 1, t1 is the time used by f1, and t2 is the time used by f2
  In Example 2, t1 is the time used by f1, and t2 is the time used by f1 and f2
  In Example 3, t1 is the time used by f1, and t2 is the time used by f1 and f3

               Example 1:      
		  ClockReset ();
		  f1();
		  t1 = ClockTotal ();
		  ClockReset ();
		  f2();
		  t2 = ClockTotal ();


	       Example 2:
		  ClockReset ();
		  f1();
		  t1 = ClockTotal ();
		  f2();
		  t2 = ClockTotal ();
		  
	       Example 3:
		  ClockReset ();
		  f1();
		  ClockStop();
		  f2();
		  ClockContinue();
		  f3()
		  t2 = ClockTotal ();

  ========================================================================= */

#ifndef _CHRONOMETER_H_INCLUDED

extern void  ClockReset ();
extern void  ClockStop ();
extern void  ClockContinue ();
extern float ClockTotal ();

#define _CHRONOMETER_H_INCLUDED
#endif
