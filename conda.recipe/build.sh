git describe --tags --dirty > $SRC_DIR/__conda_version__.txt
$PYTHON $RECIPE_DIR/format_version.py $SRC_DIR/__conda_version__.txt

rm -rf build

$PYTHON setup.py install --old-and-unmanageable
