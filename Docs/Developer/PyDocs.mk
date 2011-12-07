### PyDocs.mk --- for python documentation

# This file has basic rules for generating python docs
# from source and relocating them to current directory.

# These variables need to be defined:
#
#    PYTH_DIR = ? [directory of python source, full path or relative to here]

SOURCES.py   = $(wildcard $(PYTH_DIR)/*.py)
SOURCES.html = $(SOURCES.py:$(PYTH_DIR)/%.py=%.html)

TARGET_INDEX = dirIndex.html
#
#  Default target is "docs"
#
docs:  $(SOURCES.html) $(TARGET_INDEX)

list-docs:  
	echo py:    $(SOURCES.py)
	echo html:  $(SOURCES.html)

#
#  Rules for building docs
#
%.html:  $(PYTH_DIR)/%.py
	-d=`pwd`; cd $(PYTH_DIR); pydoc -w ./$*.py ; mv $@ $$d/$@

$(TARGET_INDEX):
	-HtmlIndex.py -f $@ $(SOURCES.html)
#
cleanHTML:
	rm -f $(SOURCES.html) $(TARGET_INDEX)
