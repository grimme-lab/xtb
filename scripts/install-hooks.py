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
import subprocess
from pathlib import Path


def main():
    # Get repository root (parent of scripts directory)
    script_dir = Path(__file__).parent.resolve()
    repo_root = script_dir.parent

    # Resolve git hooks directory using git rev-parse to support worktrees
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--git-path", "hooks"],
            cwd=repo_root,
            capture_output=True,
            text=True,
            check=True
        )
        git_hooks_path = result.stdout.strip()
        git_hooks_dir = repo_root / git_hooks_path
    except subprocess.CalledProcessError:
        print("❌ Error: Not in a git repository.", file=sys.stderr)
        print("   Make sure you're running this from the xtb repository root.", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print("❌ Error: git command not found.", file=sys.stderr)
        sys.exit(1)

    if not git_hooks_dir.is_dir():
        print(f"❌ Error: Git hooks directory not found at {git_hooks_dir}", file=sys.stderr)
        print("   Make sure you're running this from the xtb repository root.", file=sys.stderr)
        sys.exit(1)

    # Source and destination paths
    source_hook = script_dir / "commit-msg"
    dest_hook = git_hooks_dir / "commit-msg"

    if not source_hook.exists():
        print(f"❌ Error: Source hook not found at {source_hook}", file=sys.stderr)
        sys.exit(1)

    # Install hook
    print("Installing git hooks...")

    try:
        # Check if hook already exists and back it up
        if dest_hook.exists():
            # Find a unique backup filename
            backup_hook = git_hooks_dir / "commit-msg.backup"
            if backup_hook.exists():
                print(f"⚠️  Backup file already exists at {backup_hook}")
                print("   Skipping backup to preserve original hook.")
            else:
                print(f"⚠️  Existing commit-msg hook found. Backing up to {backup_hook}")
                shutil.copy2(dest_hook, backup_hook)

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
