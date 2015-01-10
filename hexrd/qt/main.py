from __future__ import absolute_import, print_function

import logging
import multiprocessing
import os
import sys
import time
import webbrowser

# use API v2 to prepare for migration to Py3/PyQt5:
import sip
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QVariant', 2)
from PyQt4 import QtCore, QtGui, uic
from PyQt4 import uic

import hexrd
from hexrd import config
from .graphicscanvas import GraphicsCanvasController
from .imageseries import get_image_series
from .materialsdialog import get_material
from .preferences import get_preferences
from .resources import resources

try:
    from IPython.qt.console.rich_ipython_widget import RichIPythonWidget
    from IPython.qt.inprocess import QtInProcessKernelManager

    # Create an in-process kernel
    kernel_manager = QtInProcessKernelManager()
    kernel_manager.start_kernel()
    kernel = kernel_manager.kernel
    kernel.gui = 'qt4'

    kernel_client = kernel_manager.client()
    kernel_client.start_channels()
except ImportError:
    pass


logger = logging.getLogger('hexrd')


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



class WhatsThisUrlLoader(QtCore.QObject):

    def eventFilter(self, target, e):
        if e.type() == QtCore.QEvent.WhatsThisClicked:
            webbrowser.open_new_tab(e.href())
            return True
        return False



class IPythonNamespaceUpdater(QtCore.QObject):

    def eventFilter(self, target, e):
        if e.type() == QtCore.QEvent.Enter:
            kernel.shell.user_ns.update(globals())
        return False



class MainController(QtGui.QMainWindow):


    @property
    def cfg(self):
        return self._cfgs[0]


    def __init__(self, log_level, cfg_file=None):
        super(MainController, self).__init__()
        # Create and display the splash screen
        splash_pix = QtGui.QPixmap(resources['hexrd.png'])
        splash = QtGui.QSplashScreen(splash_pix)#, QtCore.Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.show()
        # sleep seems necessary on linux to render the image
        time.sleep(0.01)
        QtGui.qApp.processEvents()

        self._preferences = get_preferences()

        uic.loadUi(resources['mainwindow.ui'], self)
        self.menuHelp.addAction(QtGui.QWhatsThis.createAction(self))
        self._configure_tool_buttons()

        # now that we have the ui, configure the logging widget
        add_handler(log_level, QLogStream(self.loggerTextEdit))

        self.gc_ctlr = GraphicsCanvasController(self)

        self.load_config(cfg_file)

        self.settings = QtCore.QSettings('hexrd', 'hexrd')
        self._restore_state()
        self._configure_dock_widgets()

        # load all ui elements before loading the event filters
        self._load_event_filters()

        # give the splash screen a little time to breathe
        time.sleep(1)

        self.show()

        splash.finish(self)


    def _load_ipython(self):
        def stop():
            kernel_client.stop_channels()
            kernel_manager.shutdown_kernel()

        control = RichIPythonWidget()
        self.ipythonDockWidget.setWidget(control)
        control.kernel_manager = kernel_manager
        control.kernel_client = kernel_client
        control.exit_requested.connect(stop)
        control.installEventFilter(IPythonNamespaceUpdater(self))


    def _load_event_filters(self):
        temp = WhatsThisUrlLoader(self)
        self.installEventFilter(temp)
        for k, v in self.__dict__.items():
            if isinstance(v, QtGui.QWidget):
                v.installEventFilter(temp)


    def _configure_dock_widgets(self):
        self.menuView.addAction(self.configDockWidget.toggleViewAction())
        self.menuView.addAction(self.docsDockWidget.toggleViewAction())
        try:
            self._load_ipython()
            self.menuView.addAction(self.ipythonDockWidget.toggleViewAction())
        except NameError:
            logger.info("IPython not installed, embedded console not available")
            self.ipythonDockWidget.setVisible(False)
        self.menuView.addAction(self.loggerDockWidget.toggleViewAction())


    def _configure_tool_buttons(self):
        self.imageSeriesToolButton.addAction(self.actionLoadImageSeries)
        self.imageSeriesToolButton.addAction(self.actionModifyImageSeries)
        self.imageSeriesToolButton.addAction(self.actionDeleteImageSeries)

        self.materialToolButton.addAction(self.actionAddMaterial)
        self.materialToolButton.addAction(self.actionModifyMaterial)
        self.materialToolButton.addAction(self.actionDeleteMaterial)


    def _restore_state(self):
        # main window state
        self.setWindowTitle('HEXRD')
        temp = self.settings.value('geometry')
        if temp is not None:
            self.restoreGeometry(temp)
        temp = self.settings.value('state')
        if temp is not None:
            self.restoreState(self.settings.value('state'))
        temp = self.settings.value('currentTool')
        if temp is not None:
            self.cfgToolBox.setCurrentIndex(int(temp))
        # docs
        self.docsWebView.load(QtCore.QUrl(hexrd.doc_url))
        # graphics state
        temp = self.settings.value('currentColorMap')
        if temp is not None:
            temp = self.colorMapComboBox.findText(temp)
            if temp > -1:
                self.colorMapComboBox.setCurrentIndex(temp)
        self.cmapReverseCheckBox.setChecked(
            self.settings.value('cmapReverse') == 'true'
            )
        temp = self.settings.value('imageMax')
        if temp is not None:
            self.maxSpinBox.setValue(int(temp))
        self.showOverCheckBox.setChecked(
            self.settings.value('showOver') == 'true'
            )
        temp = self.settings.value('imageMin')
        if temp is not None:
            self.minSpinBox.setValue(int(temp))
        self.showUnderCheckBox.setChecked(
            self.settings.value('showUnder') == 'true'
            )
        # view mode
        self.singleImRadioButton.setChecked(
            self.settings.value('viewSingleImage') == 'true'
            )
        self.maxImRadioButton.setChecked(
            self.settings.value('viewMaxImage') == 'true'
            )
        self.aveImRadioButton.setChecked(
            self.settings.value('viewAveImage') == 'true'
            )
        self.medianImRadioButton.setChecked(
            self.settings.value('viewMedianImage') == 'true'
            )
        self.minImRadioButton.setChecked(
            self.settings.value('viewMinImage') == 'true'
            )
        # materials
        self.showRingsCheckBox.setChecked(
             self.settings.value('showRings') == 'true'
            )
        self.showRangesCheckBox.setChecked(
             self.settings.value('showRanges') == 'true'
            )
        temp = self.settings.value('tthRanges')
        if temp is not None:
            self.tthRangesSpinBox.setValue(float(temp))
        # calibration
        temp = self.settings.value('detxStep')
        if temp is not None:
            self.detxStepSpinBox.setValue(float(temp))
        temp = self.settings.value('detyStep')
        if temp is not None:
            self.detyStepSpinBox.setValue(float(temp))
        temp = self.settings.value('detxRotStep')
        if temp is not None:
            self.detxRotStepSpinBox.setValue(float(temp))
        temp = self.settings.value('detyRotStep')
        if temp is not None:
            self.detyRotStepSpinBox.setValue(float(temp))
        temp = self.settings.value('detzRotStep')
        if temp is not None:
            self.detzRotStepSpinBox.setValue(float(temp))
        temp = self.settings.value('detDistanceStep')
        if temp is not None:
            self.detDistanceStepSpinBox.setValue(float(temp))
        temp = self.settings.value('p0DistortionStep')
        if temp is not None:
            self.p0DistortionStepSpinBox.setValue(float(temp))
        temp = self.settings.value('p1DistortionStep')
        if temp is not None:
            self.p1DistortionStepSpinBox.setValue(float(temp))
        temp = self.settings.value('p2DistortionStep')
        if temp is not None:
            self.p2DistortionStepSpinBox.setValue(float(temp))
        temp = self.settings.value('chiTiltStep')
        if temp is not None:
            self.chiTiltStepSpinBox.setValue(float(temp))


    def _save_state(self):
        # geometry
        self.settings.setValue('geometry', self.saveGeometry())
        self.settings.setValue('state', self.saveState())
        self.settings.setValue('currentTool', self.cfgToolBox.currentIndex())
        # graphics window
        self.settings.setValue(
            'currentColorMap', self.colorMapComboBox.currentText()
            )
        self.settings.setValue(
            'cmapReverse', self.cmapReverseCheckBox.isChecked()
            )
        self.settings.setValue('imageMax', self.maxSpinBox.value())
        self.settings.setValue('imageMin', self.minSpinBox.value())
        self.settings.setValue('showOver', self.showOverCheckBox.isChecked())
        self.settings.setValue('showUnder', self.showUnderCheckBox.isChecked())
        # image series
        self.settings.setValue(
            'viewSingleImage', self.singleImRadioButton.isChecked()
            )
        self.settings.setValue(
            'viewMaxImage', self.maxImRadioButton.isChecked()
            )
        self.settings.setValue(
            'viewAveImage', self.aveImRadioButton.isChecked()
            )
        self.settings.setValue(
            'viewMedianImage', self.medianImRadioButton.isChecked()
            )
        self.settings.setValue(
            'viewMinImage', self.minImRadioButton.isChecked()
            )
        # materials
        self.settings.setValue('showRings', self.showRingsCheckBox.isChecked())
        self.settings.setValue(
            'showRanges', self.showRangesCheckBox.isChecked()
            )
        self.settings.setValue('tthRanges', self.tthRangesSpinBox.value())
        # instrument calibration
        self.settings.setValue('detxStep', self.detxStepSpinBox.value())
        self.settings.setValue('detyStep', self.detyStepSpinBox.value())
        self.settings.setValue('detxRotStep', self.detxRotStepSpinBox.value())
        self.settings.setValue('detyRotStep', self.detyRotStepSpinBox.value())
        self.settings.setValue('detzRotStep', self.detzRotStepSpinBox.value())
        self.settings.setValue(
            'detDistanceStep', self.detDistanceStepSpinBox.value()
            )
        self.settings.setValue(
            'p0DistortionStep', self.p0DistortionStepSpinBox.value()
            )
        self.settings.setValue(
            'p1DistortionStep', self.p1DistortionStepSpinBox.value()
            )
        self.settings.setValue(
            'p2DistortionStep', self.p2DistortionStepSpinBox.value()
            )
        self.settings.setValue('chiTiltStep', self.chiTiltStepSpinBox.value())


    @QtCore.pyqtSlot()
    def on_actionDocumentation_triggered(self):
        if self._preferences.docs_view == 'dockWidget':
            self.docsDockWidget.setVisible(True)
        elif self._preferences.docs_view == 'webbrowser':
            webbrowser.open_new_tab(hexrd.doc_url)


    @QtCore.pyqtSlot()
    def on_actionAbout_triggered(self):
        QtGui.QMessageBox.about(
            self, 'About HEXRD',
"HEXRD provides a collection of resources for analysis of x-ray diffraction "
"data, especially high-energy x-ray diffraction. HEXRD is comprised of a "
"library and API for writing scripts, a command line interface, and an "
"interactive graphical user interface.\n\n"
"HEXRD is an open-source project originally conceived by Joel Bernier, and "
"developed by Joel Bernier, Darren Dale, and Donald Boyce, et.al."
            )


    @QtCore.pyqtSlot()
    def on_actionAddMaterial_triggered(self):
        get_material(ui=True)


    @QtCore.pyqtSlot()
    def on_actionLoadConfiguration_triggered(self):
        if self.cfg.dirty:
            QMessageBox = QtGui.QMessageBox
            confirm = QMessageBox()
            confirm.setText('Configuration has been modified.')
            confirm.setInformativeText(
"""Do you want to save your changes before loading a new configuration?"""
                )
            confirm.setStandardButtons(
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel
                )
            confirm.button(QMessageBox.Discard).setText("Discard")
            confirm.setDefaultButton(QMessageBox.Cancel)
            ret = confirm.exec_()
            if ret == QMessageBox.Save:
                if self.save_config() is False:
                    return
            elif ret == QMessageBox.Cancel:
                return
        temp = QtGui.QFileDialog.getOpenFileName(
            self, 'Load Configuration', self.cfg.working_dir,
            'YAML files (*.yml)'
            )
        if temp:
            self.load_config(temp)


    @QtCore.pyqtSlot()
    def on_actionLoadImageSeries_triggered(self):
        get_image_series(self.cfg, self.cfg.analysis_name, ui=True)


    @QtCore.pyqtSlot()
    def on_actionPreferences_triggered(self):
        self._preferences = get_preferences(ui=True)


    def on_analysisNameLineEdit_editingFinished(self):
        self.cfg.analysis_name = str(self.analysisNameLineEdit.text())


    @QtCore.pyqtSlot()
    def on_changeWorkingDirButton_clicked(self):
        temp = QtGui.QFileDialog.getExistingDirectory(
            parent=self, caption='booya', directory=self.cfg.working_dir
            )
        if temp:
            self.workingDirLineEdit.setText(temp)
            self.cfg.working_dir = temp


    @QtCore.pyqtSlot(float)
    def on_energySpinBox_valueChanged(self, val):
        try:
            self.wavelengthSpinBox.blockSignals(True)
            self.wavelengthSpinBox.setValue(12.39842/val)
        finally:
            self.wavelengthSpinBox.blockSignals(False)


    @QtCore.pyqtSlot(int)
    def on_multiprocessingSpinBox_valueChanged(self, val):
        if self.cfg.multiprocessing != val:
            self.cfg.multiprocessing = val


    @QtCore.pyqtSlot(float)
    def on_wavelengthSpinBox_valueChanged(self, val):
        try:
            self.energySpinBox.blockSignals(True)
            self.energySpinBox.setValue(12.39842/val)
        finally:
            self.energySpinBox.blockSignals(False)


    def closeEvent(self, event):
        if self._cfg_file is not None and self.cfg.dirty:
            QMessageBox = QtGui.QMessageBox
            confirm = QMessageBox()
            confirm.setText('Configuration has been modified.')
            confirm.setInformativeText('Do you want to save your changes?')
            confirm.setStandardButtons(
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel
                )
            confirm.setDefaultButton(QMessageBox.Cancel)
            ret = confirm.exec_()
            if ret == QMessageBox.Save:
                if self.save_config() is False:
                    event.ignore()
                    return
            elif ret == QMessageBox.Cancel:
                event.ignore()
                return
            else:
                self._save_state()
                event.accept()


    def load_config(self, filename):
        self._cfg_file = filename
        self._cfgs = config.open(filename)

        # general
        self.analysisNameLineEdit.setText(self.cfg.analysis_name)
        self.workingDirLineEdit.setText(self.cfg.working_dir)
        self.multiprocessingSpinBox.setMaximum(multiprocessing.cpu_count())
        self.multiprocessingSpinBox.setValue(self.cfg.multiprocessing)


    @QtCore.pyqtSlot(name='on_actionSaveConfiguration_triggered')
    def save_config(self):
        temp = QtGui.QFileDialog.getSaveFileName(
            parent=self, caption='Save Configuration',
            directory=self.cfg.working_dir, filter='YAML files (*.yml)'
            )
        if temp:
            config.save(self._cfgs, temp)
            return True
        return False


    @QtCore.pyqtSlot(name='on_actionLoadCalibration_triggered')
    @QtCore.pyqtSlot(name='on_actionLoadMaterials_triggered')
    @QtCore.pyqtSlot(name='on_actionModifyImageSeries_triggered')
    @QtCore.pyqtSlot(name='on_actionDeleteImageSeries_triggered')
    @QtCore.pyqtSlot(name='on_actionModifyMaterial_triggered')
    @QtCore.pyqtSlot(name='on_actionDeleteMaterial_triggered')
    @QtCore.pyqtSlot(name='on_actionSaveCalibration_triggered')
    @QtCore.pyqtSlot(name='on_actionSaveMaterials_triggered')
    @QtCore.pyqtSlot(name='on_actionPowderBinnedFit_triggered')
    @QtCore.pyqtSlot(name='on_actionPowderDirectFit_triggered')
    @QtCore.pyqtSlot(name='on_actionSingleCrystalFit_triggered')
    @QtCore.pyqtSlot(name='on_actionCake_triggered')
    @QtCore.pyqtSlot(name='on_actionDarkImage_triggered')
    @QtCore.pyqtSlot(name='on_actionPolarRebin_triggered')
    @QtCore.pyqtSlot(name='on_actionFindOrientations_triggered')
    @QtCore.pyqtSlot(name='on_actionFitGrains_triggered')
    def not_implemented(self):
        dlg = QtGui.QMessageBox.information(
            self, 'Not Implemented',
"""The requested feature has not been implemented.

Please consider filing an issue report at
http://github.com/praxes/hexrd/issues, if one does not already exist.""",
            buttons=QtGui.QMessageBox.Ok
            )



def execute(args):
    app = QtGui.QApplication(sys.argv)
    app.setApplicationName('HEXRD')

    # configure logging
    if args.debug:
        log_level = logging.DEBUG
        add_handler(log_level)
    else:
        log_level = logging.CRITICAL if args.quiet else logging.INFO

    ctlr = MainController(log_level, args.config)

    sys.exit(app.exec_())
