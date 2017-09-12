# RTA library

#### configure

In the `include` folder, `rta.h` includes a configuration header file
to eventually override its default definitions.
This file must be created by users of the library, and named **`rta_configuration.h`**.
The files `rta_configuration_basic.h` and `rta_configuration_console.h` located
in the `include` folder are examples of what it should look like.

#### build the doc

`~/rta-lib$ cd docs`
`~/rta-lib/docs$ doxygen Doxyfile`

then open in your favourite browser

`~/rta-lib/docs$ cd html`
`~/rta-lib/docs/html$ open index.html`

#### license

All the files in RTA are released under the BSD 3-clause license, except
`statistics/rta_cca.h` and `statistics/rta_cca.c` which depend on the GSL
(GNU Scientific Library), and are therefore released under the same GPL v3
license.