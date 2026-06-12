#!/bin/bash

# Script to create conda/mamba environments from .yaml files in the current directory
# Usage: ./setup-envs.sh [--force]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORCE=false

# Parse command line arguments
if [[ "$1" == "--force" ]]; then
    FORCE=true
fi

# Determine which package manager to use
if command -v mamba &> /dev/null; then
    PKG_MANAGER="mamba"
else
    PKG_MANAGER="conda"
fi

echo "Using $PKG_MANAGER as package manager"
echo "Script directory: $SCRIPT_DIR"
echo "Force mode: $FORCE"
echo ""

# Find all .yaml files in the script directory
yaml_files=("$SCRIPT_DIR"/*.yaml)

if [[ ! -e "${yaml_files[0]}" ]]; then
    echo "No .yaml files found in $SCRIPT_DIR"
    exit 1
fi

# Process each .yaml file
for yaml_file in "${yaml_files[@]}"; do
    if [[ ! -f "$yaml_file" ]]; then
        continue
    fi
    
    # Extract environment name from the yaml file
    env_name=$(grep "^name:" "$yaml_file" | head -1 | awk '{print $2}')
    
    if [[ -z "$env_name" ]]; then
        echo "Warning: Could not extract environment name from $yaml_file, skipping..."
        continue
    fi
    
    echo "Processing environment: $env_name"
    
    # Check if environment exists
    if $PKG_MANAGER env list | grep -q "^$env_name "; then
        if [[ "$FORCE" == true ]]; then
            echo "  Removing existing environment: $env_name"
            $PKG_MANAGER env remove -n "$env_name" -y
            echo "  Creating environment from $yaml_file"
            $PKG_MANAGER env create -f "$yaml_file" -y
            echo "  ✓ Environment $env_name created successfully"
        else
            echo "  ✓ Environment $env_name already exists, skipping"
        fi
    else
        echo "  Creating environment from $yaml_file"
        $PKG_MANAGER env create -f "$yaml_file" -y
        echo "  ✓ Environment $env_name created successfully"
    fi
    
    echo ""
done

echo "Done! All environments have been processed."
