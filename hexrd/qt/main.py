import logging
import os
import sys
import time

# use API v2 to prepare for migration to Py3/PyQt5:
import sip
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QVariant', 2)

from PyQt4.QtCore import Qt, QObject
from PyQt4.QtGui import qApp, QApplication, QMainWindow, QPixmap, QSplashScreen
from PyQt4.uic import loadUi

from hexrd.qt.resources import image_files, ui_files


class QLogStream(object):

    def __init__(self, dest):
        self._dest = dest

    def write(self, val):
        # append adds its own line endings, need to rstrip
        self._dest.append(val.rstrip())



def add_handler(log_level, stream=None):
    logger = logging.getLogger('hexrd')
    h = logging.StreamHandler(stream)
    h.setLevel(log_level)
    h.setFormatter(
        logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
        )
    logger.addHandler(h)


class MainController(QObject):


    logger = logging.getLogger('hexrd')


    def __init__(self, log_level):
        # Create and display the splash screen
        splash_pix = QPixmap(image_files['splash'])
        splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.show()
        qApp.processEvents()

        # give the splash screen a little time to breathe
        time.sleep(2)

        self.ui = ui = loadUi(ui_files['main_window'])
        ui.setWindowTitle('HEXRD')

        # now that we have the ui, configure the logging widget
        add_handler(log_level, QLogStream(ui.loggerTextEdit))

        ui.show()
        splash.finish(ui)



def execute(args):
    app = QApplication(sys.argv)
    app.setApplicationName('HEXRD')

    # configure logging
    if args.debug:
        log_level = logging.DEBUG
        add_handler(log_level)
    else:
        log_level = logging.CRITICAL if args.quiet else logging.INFO

    ctlr = MainController(log_level)

    sys.exit(app.exec_())
