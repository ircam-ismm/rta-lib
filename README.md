# RTA library

#### build the doc

`~/rta$ cd docs`
`~/rta/docs$ doxygen Doxyfile`

then open in your favourite browser

`~/rta/docs$ cd html`
`~/rta/docs/html$ open index.html`

#### configure

In the `include` folder, `rta.h` includes a configuration header file
to override its default values.  
The default configuration file is `rta_configuration.h`,
but `rta_configuration_console.h` could be included instead.  
