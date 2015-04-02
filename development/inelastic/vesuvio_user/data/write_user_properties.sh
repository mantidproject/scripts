#!/bin/bash
#
# Write out a user properties file for testing
#
PROPS_FILE=$HOME/.mantid/Mantid.user.properties

data_dir=$(dirname $(readlink -f $0))
algorithm_dir=$(dirname $data_dir)/vesuvio/algorithms
# Backup current file
if [ -f $PROPS_FILE ]; then
  mv -f $PROPS_FILE $PROPS_FILE.bak
fi

# Write out test file
echo "UpdateInstrumentDefinitions.OnStartup=0" > $PROPS_FILE
echo "usagereports.enabled=0" >> $PROPS_FILE
echo "user.python.plugins.directories = $algorithm_dir" >> $PROPS_FILE
echo "datasearch.directories = $data_dir" >> $PROPS_FILE
