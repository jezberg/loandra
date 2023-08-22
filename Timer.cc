/*!
 * \author Jeremias Berg - jeremiasberg@hmail.com
 *
 * @section LICENSE
 *   Loandra, Copyright (c) 2018, Jeremias Berg, Emir Demirovic, Peter Stuckey
 *  Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "Timer.h"

using namespace openwbo;

time_t Timer::timeSinceStart() {
  return time(NULL) - time_start;
}

void Timer::start_timer() {
  hasStarted = true;
  time_start = time(NULL);
}

bool Timer::terminate () {
  if (!hasStarted || lim_cg_sec < 0) {
    return false;
  }
  time_t cur = time(NULL);
  return (cur - time_start) > lim_cg_sec;
}


 









