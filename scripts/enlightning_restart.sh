#!/bin/bash

mpirun -np $1 ./main input.txt enlightning.restart $2
