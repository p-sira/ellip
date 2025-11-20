#!/bin/bash

set -e
trap "echo 'Interrupted'; exit 130" INT

# Generate error report
cargo run --example generate_error_report

# Inject README
cargo run --example generate_test_summary

# Generate figures
mkdir figures/bin
BINS=($(find ellip-plot-graph/src/bin -maxdepth 1 -type f -name "*.rs" -exec basename {} .rs \;))

MAX_ATTEMPTS=3

for BIN in "${BINS[@]}"; do
    echo "Running bin: $BIN"

    ATTEMPTS=0
    until [ $ATTEMPTS -ge $MAX_ATTEMPTS ]; do
        if timeout 180s cargo run -p ellip-plot-graph --bin "$BIN"; then
            echo "Plot succeeded"
            break
        else
            ATTEMPTS=$((ATTEMPTS+1))
            echo "Attempt $ATTEMPTS failed. Retrying..."
        fi
    done

    if [ $ATTEMPTS -eq $MAX_ATTEMPTS ]; then
        echo "Failed to generate $BIN — exiting."
        exit 1
    fi
done

git add .
git commit -m "CI, DOC: Generate figures and test reports."
