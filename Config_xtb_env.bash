#!/bin/bash
# run this script to set up a xtb environment
# requirements: $XTBHOME is set to `pwd`
if [ -z "${XTBHOME}" ]; then
   XTBHOME="$(cd -P "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# set up path for xtb, using the xtb directory and the users home directory
XTBPATH=${XTBHOME}:${HOME}

# to include the documentation we include our man pages in the users manpath
MANPATH=${MANPATH}:${XTBHOME}/man

# finally we have to make the binaries and scripts accessable
PATH=${PATH}:${XTBHOME}/bin:${XTBHOME}/scripts

export PATH XTBPATH MANPATH
