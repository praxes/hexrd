import logging
import os
import sys
import time

# use API v2 to prepare for migration to Py3/PyQt5:
import sip
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QVariant', 2)

from PyQt4.QtCore import Qt
from PyQt4.QtGui import QApplication, QMainWindow, QPixmap, QSplashScreen
from PyQt4.uic import loadUi

from hexrd.qt.resources import image_files, ui_files


class QLogStream(object):

    def __init__(self, dest):
        self._dest = dest

    def write(self, val):
        # append adds its own line endings, need to rstrip
        self._dest.append(val.rstrip())


def execute(args):
    app = QApplication(sys.argv)
    app.setApplicationName('HEXRD')

    # Create and display the splash screen
    splash_pix = QPixmap(image_files['splash'])
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    app.processEvents()

    # give the splash screen a little time to breathe
    time.sleep(2)

    ui = loadUi(ui_files['main_window'])
    ui.setWindowTitle('HEXRD')

    # configure logging to the log window
    # TODO: clean up
    qstream = QLogStream(ui.loggerTextEdit)
    log_level = logging.DEBUG if args.debug else logging.INFO
    if args.quiet:
        log_level = logging.ERROR
    logger = logging.getLogger('hexrd')
    logger.setLevel(log_level)
    qh = logging.StreamHandler(qstream)
    qh.setLevel(logging.CRITICAL if args.quiet else log_level)
    qh.setFormatter(
        logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
        )
    logger.addHandler(qh)
    if args.debug:
        # also log to the console, just in case
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(
            logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
            )
        logger.addHandler(ch)

    ui.show()
    splash.finish(ui)
    logger.debug("hexrd loaded")
    sys.exit(app.exec_())
