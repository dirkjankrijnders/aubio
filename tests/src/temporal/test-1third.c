#include <aubio.h>

int main(){

  uint_t samplerate = 48000;
     f = new_aubio_filter_1thirdoctave (samplerate, 13500);

          
    del_aubio_filter (f);

    
  }

  return 0;
}

