#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Specify a config file"
    exit 1
fi

datander.sh $1
datruf.sh $1
dacmaster.sh $1
dalayout.sh $1
