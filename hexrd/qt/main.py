import logging
import multiprocessing
import os
import sys
import time

# use API v2 to prepare for migration to Py3/PyQt5:
import sip
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QVariant', 2)

from PyQt4.QtCore import Qt, QObject, QSettings
from PyQt4.QtGui import (
    qApp, QApplication, QFileDialog, QMainWindow, QPixmap, QSplashScreen
    )
from PyQt4.uic import loadUi

from matplotlib import cm

from hexrd import config
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


class GraphicsCanvasController(QObject):

    def __init__(self, ui):
        self.ui = ui

        cmaps = sorted(i[:-2] for i in dir(cm) if i.endswith('_r'))
        ui.colorMapComboBox.addItems(cmaps)
        ui.colorMapComboBox.setCurrentIndex(cmaps.index('gist_heat'))

        ui.minSpinBox.valueChanged[int].connect(
            ui.maxSpinBox.setMinimum
            )
        ui.maxSpinBox.valueChanged[int].connect(
            ui.minSpinBox.setMaximum
            )


class MainController(QMainWindow):


    def __init__(self, log_level, cfg=None):
        super(MainController, self).__init__()
        # Create and display the splash screen
        splash_pix = QPixmap(image_files['splash'])
        splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.show()
        qApp.processEvents()

        # give the splash screen a little time to breathe
        time.sleep(2)

        self.ui = ui = loadUi(ui_files['main_window'], self)
        ui.setWindowTitle('HEXRD')

        # now that we have the ui, configure the logging widget
        add_handler(log_level, QLogStream(ui.loggerTextEdit))

        self.load_cfg(cfg)

        self.gc_ctlr = GraphicsCanvasController(ui)

        ui.actionExit.triggered.connect(self.close)
        ui.changeWorkingDirButton.clicked.connect(self.change_working_dir)
        ui.multiprocessingSpinBox.valueChanged[int].connect(
            lambda val: setattr(self.cfg, 'multiprocessing', val)
            )

        self.settings = QSettings('hexrd', 'hexrd')
        try:
            self.restoreGeometry(self.settings.value('geometry'))
            self.restoreState(self.settings.value('state'))
            self.cfgToolBox.setCurrentIndex(
                self.settings.value('currentTool', type=int)
                )
        except TypeError:
            raise
        ui.show()
        splash.finish(ui)


    def closeEvent(self, event):
        self.settings.setValue('geometry', self.saveGeometry())
        self.settings.setValue('state', self.saveState())
        self.settings.setValue('currentTool', self.cfgToolBox.currentIndex())
        event.accept()


    def load_cfg(self, cfg):
        self.cfg = cfg
        ui = self.ui

        # general
        ui.analysisNameLineEdit.setText(self.cfg.analysis_name)

        ui.workingDirLineEdit.setText(self.cfg.working_dir)

        ui.multiprocessingSpinBox.setMaximum(multiprocessing.cpu_count())
        ui.multiprocessingSpinBox.setValue(self.cfg.multiprocessing)


    def change_working_dir(self):
        temp = QFileDialog.getExistingDirectory(
            parent=self.ui, caption='booya', directory=self.cfg.working_dir
            )
        if temp:
            self.ui.workingDirLineEdit.setText(temp)
            self.cfg.working_dir = temp
            assert self.cfg.working_dir == temp



def execute(args):
    app = QApplication(sys.argv)
    app.setApplicationName('HEXRD')

    # configure logging
    if args.debug:
        log_level = logging.DEBUG
        add_handler(log_level)
    else:
        log_level = logging.CRITICAL if args.quiet else logging.INFO

    cfg = config.open(args.config)[0]

    ctlr = MainController(log_level, cfg)

    sys.exit(app.exec_())
