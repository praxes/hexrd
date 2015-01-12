from PyQt4 import QtCore


class WhatsThisUrlLoader(QtCore.QObject):
    def eventFilter(self, target, e):
        if e.type() == QtCore.QEvent.WhatsThisClicked:
            webbrowser.open_new_tab(e.href())
            return True
        return False
