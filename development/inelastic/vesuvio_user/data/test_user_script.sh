#!/bin/bash

test_dir=$(dirname $(readlink -f $0))
script_dir=$(dirname $test_dir)

# Write user properties
$($test_dir/write_user_properties.sh)

export PYTHONPATH=$script_dir:$PYTHONPATH
if [ -z "$MANTIDPATH" ]; then
  export MANTIDPATH=/opt/Mantid/bin
  export PYTHONPATH=/opt/Mantid/bin:$PYTHONPATH
  echo "Python path set to $PYTHONPATH"
fi

echo "Running user script..."
echo
python $script_dir/vesuvio_user_script.py
