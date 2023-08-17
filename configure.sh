#!/bin/sh
	echo "Making Cadical" 
	(cd cadical && ./configure && make)
    echo "MaxPre2"
	(cd maxpre2 && make lib)