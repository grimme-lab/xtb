#!/usr/bin/env python3
"""
Install git hooks for xtb development

This file is part of xtb.
SPDX-Identifier: LGPL-3.0-or-later
"""

import os
import sys
import shutil
import stat
from pathlib import Path


def main():
    # Get repository root (parent of scripts directory)
    script_dir = Path(__file__).parent.resolve()
    repo_root = script_dir.parent

    # Check if we're in a git repository
    git_hooks_dir = repo_root / ".git" / "hooks"
    if not git_hooks_dir.is_dir():
        print("❌ Error: .git/hooks directory not found.")
        print("   Make sure you're running this from the xtb repository root.")
        sys.exit(1)

    # Source and destination paths
    source_hook = script_dir / "commit-msg"
    dest_hook = git_hooks_dir / "commit-msg"

    if not source_hook.exists():
        print(f"❌ Error: Source hook not found at {source_hook}")
        sys.exit(1)

    # Install hook
    print("Installing git hooks...")

    try:
        shutil.copy2(source_hook, dest_hook)

        # Make executable on Unix-like systems
        if os.name != "nt":  # Not Windows
            st = os.stat(dest_hook)
            os.chmod(dest_hook, st.st_mode | stat.S_IEXEC | stat.S_IXGRP | stat.S_IXOTH)

        print("✅ Git hooks installed successfully!")
        print("")
        print("📝 All commits must now be signed off using:")
        print('   git commit -s -m "your message"')
        print("")
        print("   Or amend existing commits:")
        print("   git commit --amend --signoff")

    except Exception as e:
        print(f"❌ Error installing hooks: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
