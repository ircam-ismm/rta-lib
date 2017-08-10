# RTA library

#### configure

In the `include` folder, `rta.h` includes a configuration header file
to override its default values.  
The default configuration file is `rta_configuration.h`,
but `rta_configuration_console.h` could be included instead.  

#### build the doc

`~/rta-lib$ cd doc`  
`~/rta-lib/doc$ doxygen Doxyfile`

then open in your favourite browser

`~/rta-lib/doc$ cd html`  
`~/rta-lib/doc/html$ open index.html`

#### license

All the files in RTA are released under the BSD 3-clause license, except
`statistics/rta_cca.h` and `statistics/rta_cca.c` which depend on the GSL
(GNU Scientific Library), and are therefore released under the same GPL v3
license.