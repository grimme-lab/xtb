#!/usr/bin/env bash
# This file is part of xtb.
# SPDX-Identifier: LGPL-3.0-or-later
#
# xtb is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# xtb is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with xtb.  If not, see <https://www.gnu.org/licenses/>.

if [ -z "$MESON_DIST_ROOT" ]; then
    echo "Error: This script must be run via meson dist."
    exit 1
fi

echo "Setting up symlinks for nested subprojects in distribution tarball..."

cd "$MESON_DIST_ROOT/subprojects" || exit 1

link_subproject() {
    local parent="$1"
    local dependency="$2"

    if [ -d "$parent" ] && [ -d "$dependency" ]; then
        echo "  -> Symlinking $dependency into $parent/subprojects"
        mkdir -p "$parent/subprojects"
        ln -sfn "../../$dependency" "$parent/subprojects/$dependency"
    fi
}

link_project_alias() {
    local source="$1"
    local alias="$2"

    if [ -d "$source" ] && { [ ! -e "$alias" ] || [ -L "$alias" ]; }; then
        echo "  -> Symlinking $source as top-level $alias subproject"
        ln -sfn "$source" "$alias"
    fi
}

# Meson uses the wrap/subproject name cpx, while CMake searches for cpcmx.
link_project_alias cpx cpcmx

# tblite is a direct xTB dependency, but in the xTB distribution tarball its
# CMake subprojects live one directory above tblite/subprojects.
link_subproject tblite toml-f
link_subproject tblite multicharge
link_subproject tblite dftd4
link_subproject tblite s-dftd3
link_subproject tblite mstore

# Same pattern as tblite: dftd4 and toml-f need these nested sources for
# CMake/offline builds.
link_subproject dftd4 mstore
link_subproject toml-f test-drive

# CPCM-X is vendored as cpx in Meson, but its CMake build has its own nested
# dependency lookup.
link_subproject cpx numsa
link_subproject cpx toml-f
link_subproject cpx test-drive
link_subproject numsa test-drive

# Optional JSON support in mctc-lib/jonquil needs these sources offline.
link_subproject mctc-lib toml-f
link_subproject mctc-lib jonquil
link_subproject jonquil test-drive

echo "Symlink setup complete."
