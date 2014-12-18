from matplotlib import cm

from PyQt4 import QtCore


class GraphicsCanvasController(QtCore.QObject):

    def __init__(self, ui):
        QtCore.QObject.__init__(self)
        self.ui = ui

        cmaps = sorted(i[:-2] for i in dir(cm) if i.endswith('_r'))
        ui.colorMapComboBox.addItems(cmaps)

        ui.minSpinBox.valueChanged[int].connect(
            ui.maxSpinBox.setMinimum
            )
        ui.maxSpinBox.valueChanged[int].connect(
            ui.minSpinBox.setMaximum
            )
