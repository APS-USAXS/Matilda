#!/bin/bash

# run the tiled server in a screen session
SESSION_NAME=matilda
N_LINES=5000

screen \
    -dm \
    -S "${SESSION_NAME}" \
    -h "${N_LINES}" \
    ./st-matilda.sh
