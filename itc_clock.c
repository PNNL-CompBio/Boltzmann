#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/time.h>

int64_t itc_clock() {
  struct timeval tv;
  struct timezone *tz;
  int64_t clock_v;
  int64_t sec_v;
  int success;
  tz = NULL;
  success = gettimeofday(&tv,tz);
  /*
  clock_v = (int64_t)tv.tv_sec * 1000000L + (int64_t)tv.tv_usec;
  */
  sec_v = (int64_t)tv.tv_sec;
  clock_v =   (sec_v << 20);
  clock_v -=  (sec_v << 16);
  clock_v +=  (sec_v << 14); 
  clock_v +=  (sec_v << 9);
  clock_v +=  (sec_v << 6);
  clock_v +=  (int64_t)tv.tv_usec;
  return(clock_v);
}
