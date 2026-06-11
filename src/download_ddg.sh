#!/bin/sh
# clone thermompnn ddG data into data/ThermoMPNN (we only use the test split)
set -e
d="$(dirname "$0")/../data/ThermoMPNN"
git clone https://github.com/Kuhlman-Lab/ThermoMPNN "$d"
