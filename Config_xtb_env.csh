#!/bin/csh
# run this script to set up a xtb environment
# requirements: $XTBHOME is set to `pwd`
if ( ! ${?XTBHOME} ) then
   setenv XTBHOME "`dirname $0`"
endif

# set up path for xtb, using the xtb directory and the users home directory
setenv XTBPATH ${XTBHOME}:${HOME}

# to include the documentation we include our man pages in the users manpath
setenv MANPATH ${MANPATH}:${XTBHOME}/man

# finally we have to make the binaries and scripts accessable
setenv PATH ${PATH}:${XTBHOME}/bin:${XTBHOME}/scripts
