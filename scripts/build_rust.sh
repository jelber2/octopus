#!/usr/bin/env bash
set -e

# Restore Rust source files from git (working tree may be empty in fresh shells)
git ls-tree -r --name-only HEAD \
  | grep -E "^(rust_src|Cargo\.toml|Cargo\.lock|replit\.md)" \
  | while IFS= read -r f; do
      dir=$(dirname "$f")
      [ "$dir" != "." ] && mkdir -p "$dir"
      git cat-file blob "HEAD:$f" > "$f"
    done

# Ensure cargo is on PATH
export PATH="$HOME/.nix-profile/bin:$PATH"

cargo build --release 2>&1
echo ""
echo "=== Build complete ==="
echo ""
./target/release/octopus --help
