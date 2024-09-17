#!/bin/sh
	echo "Configuring Cadical" 
	(cd cadical && ./configure)
    echo "MaxPre2"
	(cd maxpre2 && make lib)
