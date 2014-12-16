from PyQt4 import QtGui, QtCore, uic

from .resources import resources


class Preferences(object):


    @property
    def docs_view(self):
        return self._docs_view


    def __init__(self):
        settings = QtCore.QSettings('hexrd', 'preferences')
        self._docs_view = settings.value('docsView', 'dockWidget')



class PreferencesController(QtGui.QDialog):


    def __init__(self):
        super(PreferencesController, self).__init__()
        uic.loadUi(resources['preferences.ui'], self)
        self.load_settings()


    def load_settings(self):
        prefs = Preferences()
        # TODO: fix this when buttongroup is available via PyQt
        self.docsDockableWindowButton.setChecked(
            prefs.docs_view == 'dockWidget'
            )
        self.docsWebbrowserButton.setChecked(
            prefs.docs_view == 'webbrowser'
            )


    @QtCore.pyqtSlot(name='on_buttonBox_accepted')
    def save_settings(self):
        settings = QtCore.QSettings('hexrd', 'preferences')
        if self.docsDockableWindowButton.isChecked():
            settings.setValue('docsView', 'dockWidget')
        elif self.docsWebbrowserButton.isChecked():
            settings.setValue('docsView', 'webbrowser')



def get_preferences(ui=False):
    if ui:
        ctlr = PreferencesController()
        ctlr.exec_()
    return Preferences()
