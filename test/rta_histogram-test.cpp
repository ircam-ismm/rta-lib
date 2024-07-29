// unit test for rta histogram using catch v1 framework (run by pipo test target)
#include "catch.hpp"
#include "rta_histogram.h"

TEST_CASE("rta_histogram")
{
  rta_histogram_params_t hist;
  rta_histogram_init(&hist);

  rta_real_t data[] = {2, 3, 3, 4, 4, 4};
  int size = sizeof(data) / sizeof(*data);
    
  hist.nhist = 3;
  rta_real_t output[hist.nhist];
  rta_real_t binout[hist.nhist];

  SECTION("Unweighted")
  {
    rta_histogram_stride(&hist, data, 1, size, output, 1, binout, 1);

    // verify return bin x coordinates
    CHECK(binout[0] == 2);
    CHECK(binout[1] == 3);
    CHECK(binout[2] == 4);

    // verify returned counts
    CHECK(output[0] == 1);
    CHECK(output[1] == 2);
    CHECK(output[2] == 3);

    CHECK(hist.lo == 2);
    CHECK(hist.hi == 4);
  }

  SECTION("Weighted")
  {
    rta_real_t weights[] = {1, 1, 1, 0, 0, 0};
    
    rta_histogram_weighted_stride(&hist, data, 1, size, weights, 1, output, 1, binout, 1);

    // verify return bin x coordinates
    CHECK(binout[0] == 2);
    CHECK(binout[1] == 3);
    CHECK(binout[2] == 4);

    // verify returned counts
    CHECK(output[0] == 1);
    CHECK(output[1] == 2);
    CHECK(output[2] == 0);

    CHECK(hist.lo == 2);
    CHECK(hist.hi == 4);
  }

  SECTION("No data")
  {
    // guard values to check nothing has been written
    rta_real_t data[1]  = {97};
    rta_real_t output[] = {98, 99, 100};
    rta_real_t binout[] = {99, 100, 101};

    rta_histogram_stride(&hist, data, 1, 0, output, 1, binout, 1);

    // check that lo/hi is not written, output and bins are cleared
    CHECK(hist.lo == 0);
    CHECK(hist.hi == 0);

    CHECK(output[0] == 0);
    CHECK(output[1] == 0);
    CHECK(output[2] == 0);

    CHECK(binout[0] == 0);
    CHECK(binout[1] == 0);
    CHECK(binout[2] == 0);
  }
}

/** EMACS **
 * Local variables:
 * mode: c++
 * c-basic-offset:2
 * End:
 */
