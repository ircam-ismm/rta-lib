# RTA library

#### configure

In the `include` folder, `rta.h` includes a configuration header file
to override its default values.  
The default configuration file is `rta_configuration.h`,
but `rta_configuration_console.h` could be included instead.  

#### build the doc

`~/rta$ cd docs`
`~/rta/docs$ doxygen Doxyfile`

then open in your favourite browser

`~/rta/docs$ cd html`
`~/rta/docs/html$ open index.html`

#### license

All the files in RTA are released under the BSD 3-clause license, except
`statistics/rta_cca.h` and `statistics/rta_cca.c` which depend on the
GNU Scientific Library (gsl), and are therefore released under the same GPL v3
license.  
`signal/psy.h`, `signal/psy.c`, `signal/rta_yin.h` and `signal/rta_yin.c` are all
BSD 3-clause, but there is a running patent on their commercial use.
