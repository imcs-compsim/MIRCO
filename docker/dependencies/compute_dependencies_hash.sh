#!/bin/bash

# Exit the script at the first failure
set -e

find dependencies/trilinos docker/dependencies -type f -exec sha1sum {} \; | sort | sha1sum | cut -c -8
