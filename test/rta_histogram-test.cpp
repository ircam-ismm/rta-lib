
#include "catch.hpp"
#include "rta_histogram.h"

TEST_CASE("rta_histogram")
{
  rta_histogram_params_t hist;
  rta_histogram_init(&hist);

  SECTION("Data")
  {
    rta_real_t data[] = {2, 3, 3, 4, 4, 4};
    int size = sizeof(data) / sizeof(*data);
    
    hist.nhist = 3;
    rta_real_t output[hist.nhist];
    rta_real_t binout[hist.nhist];

    rta_histogram_stride(&hist, data, 1, size, output, 1, binout, 1);

    // verify return bin x coordinates
    CHECK((float) binout[0] == 2);
    CHECK((float) binout[1] == 3);
    CHECK((float) binout[2] == 4);

    // verify returned counts
    CHECK((float) output[0] == 1);
    CHECK((float) output[1] == 2);
    CHECK((float) output[2] == 3);
  }
}

/** EMACS **
 * Local variables:
 * mode: c++
 * c-basic-offset:2
 * End:
 */
