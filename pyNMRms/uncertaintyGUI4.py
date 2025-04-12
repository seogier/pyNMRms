# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'designer\uncertaintyGUI.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_Uncertainty(object):
    def setupUi(self, Uncertainty):
        Uncertainty.setObjectName(_fromUtf8("Uncertainty"))
        Uncertainty.resize(1186, 729)
        self.centralwidget = QtGui.QWidget(Uncertainty)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.pbCalculateUncertainty = QtGui.QPushButton(self.centralwidget)
        self.pbCalculateUncertainty.setGeometry(QtCore.QRect(40, 440, 401, 57))
        self.pbCalculateUncertainty.setObjectName(_fromUtf8("pbCalculateUncertainty"))
        self.layoutWidget = QtGui.QWidget(self.centralwidget)
        self.layoutWidget.setGeometry(QtCore.QRect(30, 70, 561, 323))
        self.layoutWidget.setObjectName(_fromUtf8("layoutWidget"))
        self.verticalLayout = QtGui.QVBoxLayout(self.layoutWidget)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.formLayout = QtGui.QFormLayout()
        self.formLayout.setObjectName(_fromUtf8("formLayout"))
        self.label = QtGui.QLabel(self.layoutWidget)
        self.label.setObjectName(_fromUtf8("label"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.LabelRole, self.label)
        self.leT1Reported = QtGui.QLineEdit(self.layoutWidget)
        self.leT1Reported.setObjectName(_fromUtf8("leT1Reported"))
        self.formLayout.setWidget(0, QtGui.QFormLayout.FieldRole, self.leT1Reported)
        self.label_2 = QtGui.QLabel(self.layoutWidget)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_2)
        self.leT2Reported = QtGui.QLineEdit(self.layoutWidget)
        self.leT2Reported.setObjectName(_fromUtf8("leT2Reported"))
        self.formLayout.setWidget(1, QtGui.QFormLayout.FieldRole, self.leT2Reported)
        self.verticalLayout.addLayout(self.formLayout)
        self.formLayout_2 = QtGui.QFormLayout()
        self.formLayout_2.setObjectName(_fromUtf8("formLayout_2"))
        self.label_3 = QtGui.QLabel(self.layoutWidget)
        self.label_3.setObjectName(_fromUtf8("label_3"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_3)
        self.leT1Uncertainty = QtGui.QLineEdit(self.layoutWidget)
        self.leT1Uncertainty.setObjectName(_fromUtf8("leT1Uncertainty"))
        self.formLayout_2.setWidget(0, QtGui.QFormLayout.FieldRole, self.leT1Uncertainty)
        self.label_4 = QtGui.QLabel(self.layoutWidget)
        self.label_4.setObjectName(_fromUtf8("label_4"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_4)
        self.leT2Uncertainty = QtGui.QLineEdit(self.layoutWidget)
        self.leT2Uncertainty.setObjectName(_fromUtf8("leT2Uncertainty"))
        self.formLayout_2.setWidget(1, QtGui.QFormLayout.FieldRole, self.leT2Uncertainty)
        self.verticalLayout.addLayout(self.formLayout_2)
        self.formLayout_3 = QtGui.QFormLayout()
        self.formLayout_3.setObjectName(_fromUtf8("formLayout_3"))
        self.label_5 = QtGui.QLabel(self.layoutWidget)
        self.label_5.setObjectName(_fromUtf8("label_5"))
        self.formLayout_3.setWidget(0, QtGui.QFormLayout.LabelRole, self.label_5)
        self.ledT1dT = QtGui.QLineEdit(self.layoutWidget)
        self.ledT1dT.setObjectName(_fromUtf8("ledT1dT"))
        self.formLayout_3.setWidget(0, QtGui.QFormLayout.FieldRole, self.ledT1dT)
        self.label_6 = QtGui.QLabel(self.layoutWidget)
        self.label_6.setObjectName(_fromUtf8("label_6"))
        self.formLayout_3.setWidget(1, QtGui.QFormLayout.LabelRole, self.label_6)
        self.ledT2dT = QtGui.QLineEdit(self.layoutWidget)
        self.ledT2dT.setObjectName(_fromUtf8("ledT2dT"))
        self.formLayout_3.setWidget(1, QtGui.QFormLayout.FieldRole, self.ledT2dT)
        self.verticalLayout.addLayout(self.formLayout_3)
        Uncertainty.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(Uncertainty)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1186, 47))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        Uncertainty.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(Uncertainty)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        Uncertainty.setStatusBar(self.statusbar)

        self.retranslateUi(Uncertainty)
        QtCore.QMetaObject.connectSlotsByName(Uncertainty)

    def retranslateUi(self, Uncertainty):
        Uncertainty.setWindowTitle(_translate("Uncertainty", "Uncertainty", None))
        self.pbCalculateUncertainty.setText(_translate("Uncertainty", "Calculate Uncertainty", None))
        self.label.setText(_translate("Uncertainty", "T1 reported (ms)", None))
        self.label_2.setText(_translate("Uncertainty", "T2 reported (ms)", None))
        self.label_3.setText(_translate("Uncertainty", "T1 uncertainty (ms)", None))
        self.label_4.setText(_translate("Uncertainty", "T2 uncertainty (ms)", None))
        self.label_5.setText(_translate("Uncertainty", "dT1/dT (%/C)", None))
        self.label_6.setText(_translate("Uncertainty", "dT2/dT (%/C)", None))

