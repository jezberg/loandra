/*!
 * \author Jeremias Berg - jeremiasberg@gmail.com
 *  This implementation is heavily based on Open-WBO, thanks to the authors! 
 * 
 * @section LICENSE
 *  Loandra, Copyright (c) 2018 Jeremias Berg
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

#ifndef Timer_h
#define Timer_h

#include "utils/System.h"
#include "cadical/src/cadical.hpp"
#include <string> 
#include <ctime>


namespace openwbo {

class Timer : public CaDiCaL::Terminator {

public:
  // NOTE: currently the encoding is not set as an input parameter.
  Timer(int limit_cg_sec) {
    hasStarted = false;
    time_start = 0;
    lim_cg_sec = limit_cg_sec;
  }

  ~Timer() {

  }

  void start_timer();
  time_t timeSinceStart();
  // 0 if not terminate
  bool terminate ();


protected:

  bool hasStarted;
  int lim_cg_sec;
  time_t time_start;
  std::string print_timeSinceStart();

};
} // namespace openwbo

#endif
