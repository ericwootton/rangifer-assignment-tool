#!/bin/bash
# ============================================================================
# Rangifer Assignment Tool — Mac Launcher
# ============================================================================
# Double-click this file to launch the app.
# Requires R to be installed (https://cran.r-project.org/).
# ============================================================================

cd "$(dirname "$0")"

# Find Rscript
if command -v Rscript &>/dev/null; then
    RSCRIPT=Rscript
elif [ -f "/usr/local/bin/Rscript" ]; then
    RSCRIPT=/usr/local/bin/Rscript
elif [ -f "/opt/homebrew/bin/Rscript" ]; then
    RSCRIPT=/opt/homebrew/bin/Rscript
elif [ -d "/Library/Frameworks/R.framework" ]; then
    RSCRIPT=/Library/Frameworks/R.framework/Resources/bin/Rscript
else
    echo ""
    echo "ERROR: R is not installed."
    echo "Please install R from https://cran.r-project.org/"
    echo ""
    echo "Press Enter to close..."
    read
    exit 1
fi

echo "Using: $RSCRIPT"
"$RSCRIPT" launch.R

echo ""
echo "App stopped. Press Enter to close..."
read
