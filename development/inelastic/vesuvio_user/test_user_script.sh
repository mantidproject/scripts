#!/bin/bash

SCRIPT_DIR=`dirname $0`

# Copy user properties
cp $SCRIPT_DIR/Mantid.user.properties $HOME/.mantid

if [ -z "$MANTIDPATH" ]; then
  export MANTIDPATH=/opt/Mantid/bin
  export PYTHONPATH=/opt/Mantid/bin:$PYTHONPATH
  echo "Python path set to $PYTHONPATH"
fi

echo "Running user script..."
echo
python $SCRIPT_DIR/vesuvio_user_script.py
