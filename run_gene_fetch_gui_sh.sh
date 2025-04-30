#!/bin/bash

echo "=============================================="
echo "         Gene Fetch GUI Launcher"
echo "=============================================="
echo ""

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is not installed or not in your PATH."
    echo "Please install Python 3.9 or higher."
    read -p "Press Enter to continue..."
    exit 1
fi

# Run the application
echo "Starting Gene Fetch GUI..."
python3 gene_fetch_gui.py
if [ $? -ne 0 ]; then
    echo ""
    echo "Error: Gene Fetch GUI failed to run properly."
    echo "Please check the logs for more information."
    read -p "Press Enter to continue..."
    exit 1
fi

echo "Gene Fetch GUI completed successfully."
exit 0
