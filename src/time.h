#ifndef TIME_H
#define TIME_H

#include <sys/time.h>

/* get the time of day from the system clock, and store it (in seconds) */
double get_wall_time (void) {
/*#if defined(_MSC_VER)
  double t;

  t = GetTickCount();
  t = 1e-3 * t;

  return t;
#else*/
  struct timeval tm;
  struct timezone tz;

  gettimeofday(&tm, &tz);
  return((double)(tm.tv_sec) + 1e-6*(double)(tm.tv_usec));
//#endif
}

#endif
