# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'designer\multi-inputDialog.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(300, 200)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setGeometry(QtCore.QRect(20, 120, 191, 32))
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setGeometry(QtCore.QRect(31, 31, 30, 16))
        self.label.setObjectName("label")
        self.dsbX = QtWidgets.QDoubleSpinBox(Dialog)
        self.dsbX.setGeometry(QtCore.QRect(67, 31, 47, 20))
        self.dsbX.setDecimals(1)
        self.dsbX.setMinimum(0.1)
        self.dsbX.setMaximum(10.0)
        self.dsbX.setSingleStep(0.1)
        self.dsbX.setProperty("value", 1.0)
        self.dsbX.setObjectName("dsbX")
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setGeometry(QtCore.QRect(10, 10, 271, 16))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(Dialog)
        self.label_3.setGeometry(QtCore.QRect(31, 57, 30, 16))
        self.label_3.setObjectName("label_3")
        self.dsbY = QtWidgets.QDoubleSpinBox(Dialog)
        self.dsbY.setGeometry(QtCore.QRect(67, 57, 47, 20))
        self.dsbY.setDecimals(1)
        self.dsbY.setMinimum(0.1)
        self.dsbY.setMaximum(10.0)
        self.dsbY.setSingleStep(0.1)
        self.dsbY.setProperty("value", 1.0)
        self.dsbY.setObjectName("dsbY")
        self.label_4 = QtWidgets.QLabel(Dialog)
        self.label_4.setGeometry(QtCore.QRect(31, 83, 30, 16))
        self.label_4.setObjectName("label_4")
        self.dsbZ = QtWidgets.QDoubleSpinBox(Dialog)
        self.dsbZ.setGeometry(QtCore.QRect(67, 83, 47, 20))
        self.dsbZ.setDecimals(1)
        self.dsbZ.setMinimum(0.1)
        self.dsbZ.setMaximum(10.0)
        self.dsbZ.setSingleStep(0.1)
        self.dsbZ.setProperty("value", 1.0)
        self.dsbZ.setObjectName("dsbZ")

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Interpolate/Scale Image"))
        self.label.setText(_translate("Dialog", "Xscale"))
        self.label_2.setText(_translate("Dialog", "Incease/decrease number of voxels by scale factors"))
        self.label_3.setText(_translate("Dialog", "Yscale"))
        self.label_4.setText(_translate("Dialog", "Zscale"))

