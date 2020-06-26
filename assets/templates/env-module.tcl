#%Module
set prefix @prefix@

module-whatis "@description@"

prepend-path XTBPATH $prefix/@datadir@
prepend-path PATH $prefix/@bindir@
prepend-path MANPATH $prefix/@mandir@
prepend-path PKG_CONFIG_PATH $prefix/@libdir@/pkgconfig

# Only allow to load one instance of @name@
conflict @name@
