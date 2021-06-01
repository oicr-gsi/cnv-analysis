#!/bin/bash
cd $1
ls -C | xargs wc -l | sort -n | sed 's!\S*/!\t!'
