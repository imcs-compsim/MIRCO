#!/bin/bash

# Exit the script at the first failure
set -e

if [ ! -f "./create-mirco-python-venv.sh" ]; then
    echo "Please run this script from the root directory of the repository."
    exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

# Path to the python virtual environment.
PYTHON_VENV="${SCRIPT_DIR}/utilities/mirco-python-venv"

# If the virtual environment already exists, delete it.
if [ -d "$PYTHON_VENV" ]; then rm -Rf $PYTHON_VENV; fi

# Path to python
PYTHON_PATH=${1:-python3}

# Setup the virtual environment and source it.
$PYTHON_PATH -m venv "${PYTHON_VENV}"
source "${PYTHON_VENV}"/bin/activate

# Install all the modules defined in requirements.txt.
pip install --upgrade pip
pip install wheel
pip install -r requirements.txt

# Install the pre-commit hooks.
pre-commit install

set +e
