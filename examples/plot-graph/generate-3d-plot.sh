#!/bin/bash

# List of 3D plot binaries
BINARIES=(
    "ellipdinc-3d"
    "ellipeinc-3d"
    "ellipf-3d"
    "ellippiinc-3d"
    "ellippi-3d"
)

# Build and run each binary, waiting for user input
for BIN in "${BINARIES[@]}"; do
    echo "Running $BIN..."
    cargo run --bin "$BIN" --features open-html
    echo
    read -p "Press Enter to continue to the next 3D plot..." dummy
    echo
done

echo "All 3D plot binaries have been run."
