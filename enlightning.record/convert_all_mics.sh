#!/bin/bash

i=1

for f in mic-?.txt
do
	./mic $1 $2 $f mic-$i-info.txt
	((i++))
done
