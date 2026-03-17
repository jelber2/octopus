#!/usr/bin/env bash
# Build script for Octopus (Rust).
# Restores any missing Rust source files from the last git commit,
# then compiles and tests the project.

set -euo pipefail

export PATH="$HOME/.nix-profile/bin:$PATH"

WORKSPACE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$WORKSPACE"

echo "[build] Checking for missing tracked files..."
git ls-tree -r --name-only HEAD \
    | grep -E "^(rust_src|Cargo)" \
    | while IFS= read -r f; do
        if [ ! -f "$f" ]; then
            echo "[build]   restoring: $f"
            mkdir -p "$(dirname "$f")"
            git cat-file blob "HEAD:$f" > "$f"
        fi
    done
echo "[build] Source files OK."

echo "[build] Running: cargo build --release"
cargo build --release

echo "[build] Running: cargo test"
cargo test

echo "[build] All done."
