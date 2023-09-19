# -*- coding: utf-8 -*-
"""
/***************************************************************************
 ARMMDialog
                                 A QGIS plugin
 ARMM
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2023-06-29
        git sha              : $Format:%H$
        copyright            : (C) 2023 by GPN_GEO
        email                : nastyashevchenkomail@mail.ru
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

import os
from datetime import datetime, date

from PyQt5.QtCore import QStringListModel
from PyQt5.QtWidgets import QMainWindow, QListView, QAbstractItemView, QWidget, QVBoxLayout, QLabel, QDialog
from qgis.PyQt import uic
from qgis.PyQt import QtWidgets

# This loads your .ui file so that PyQt can populate your plugin with the elements from Qt Designer
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'ARMM_dialog_base.ui'))

from PyQt5 import QtCore, QtGui
from qgsfilewidget import QgsFileWidget



class Ui_MainWindow(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        # Required by Qt4 to initialize the UI
        self.setupUi(self)

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.setEnabled(True)
        MainWindow.resize(900, 900)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.centralwidget)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.label.setObjectName("label")
        self.verticalLayout_3.addWidget(self.label)
        self.cmbLic = QtWidgets.QComboBox(self.centralwidget)
        self.cmbLic.setObjectName("comboBox")
        self.verticalLayout_3.addWidget(self.cmbLic)
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setObjectName("pushButton")
        self.verticalLayout_3.addWidget(self.pushButton)
        self.verticalLayout_4.addLayout(self.verticalLayout_3)
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setEnabled(True)
        self.tabWidget.setObjectName("tabWidget")
        self.Wellpads = QtWidgets.QWidget()
        self.Wellpads.setObjectName("tab_5")
        self.pushButton_2 = QtWidgets.QPushButton(self.Wellpads)
        self.pushButton_2.setGeometry(QtCore.QRect(0, 490, 841, 28))
        self.pushButton_2.setObjectName("pushButton_2")
        self.lineEdit = QtWidgets.QLineEdit(self.Wellpads)
        self.lineEdit.setGeometry(QtCore.QRect(0, 450, 671, 22))
        self.lineEdit.setObjectName("lineEdit")
        self.listView_wp = QtWidgets.QListView(self.Wellpads)
        self.listView_wp.setGeometry(QtCore.QRect(10, 0, 311, 441))
        self.listView_wp.setObjectName("listView_3")
        self.tableView_2 = QtWidgets.QTableView(self.Wellpads)
        self.tableView_2.setGeometry(QtCore.QRect(340, 0, 501, 441))
        self.tableView_2.setObjectName("tableView_2")
        self.checkBox = QtWidgets.QCheckBox(self.Wellpads)
        self.checkBox.setGeometry(QtCore.QRect(700, 450, 121, 20))
        self.checkBox.setObjectName("checkBox")
        self.tabWidget.addTab(self.Wellpads, "")
        self.Rigs = QtWidgets.QWidget()
        self.Rigs.setObjectName("tab_4")
        self.listView_rig = QtWidgets.QListView(self.Rigs)
        self.listView_rig.setGeometry(QtCore.QRect(10, 0, 311, 441))
        self.listView_rig.setObjectName("listView_2")
        self.pushButton_3 = QtWidgets.QPushButton(self.Rigs)
        self.pushButton_3.setGeometry(QtCore.QRect(0, 490, 841, 28))
        self.pushButton_3.setObjectName("pushButton_3")
        self.lineEdit_2 = QtWidgets.QLineEdit(self.Rigs)
        self.lineEdit_2.setGeometry(QtCore.QRect(0, 450, 671, 22))
        self.lineEdit_2.setObjectName("lineEdit_2")

        self.tableView = QtWidgets.QTableView(self.Rigs)
        self.tableView.setGeometry(QtCore.QRect(340, 0, 501, 441))
        self.tableView.setObjectName("tableView")
        self.checkBox_2 = QtWidgets.QCheckBox(self.Rigs)
        self.checkBox_2.setGeometry(QtCore.QRect(700, 450, 121, 20))
        self.checkBox_2.setObjectName("checkBox_2")
        self.tabWidget.addTab(self.Rigs, "")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.tabWidget.addTab(self.tab, "")
        self.tab_3 = QtWidgets.QWidget()
        self.tab_3.setObjectName("tab_3")
        self.tabWidget.addTab(self.tab_3, "")

        self.tab_drill = QtWidgets.QWidget()
        self.tab_drill.setObjectName("tab_3")
        self.tabWidget.addTab(self.tab_drill, "")

        self.label_11 = QtWidgets.QLabel(self.tab_drill)
        self.label_11.setGeometry(QtCore.QRect(60, 30, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.tableWidget_7 = QtWidgets.QTableWidget(self.tab_drill)
        self.tableWidget_7.setGeometry(QtCore.QRect(20, 60, 161, 561))
        self.tableWidget_7.setObjectName("tableWidget_7")
        self.tableWidget_7.setColumnCount(1)
        self.tableWidget_7.setRowCount(100)
        # self.tableWidget_7.isEnabled = True
        # self.tableWidget_7.setEditTriggers(QAbstractItemView::NoEditTriggers)
        self.label_12 = QtWidgets.QLabel(self.tab_drill)
        self.label_12.setGeometry(QtCore.QRect(260, 30, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.tableWidget_8 = QtWidgets.QTableWidget(self.tab_drill)
        self.tableWidget_8.setGeometry(QtCore.QRect(210, 60, 161, 561))
        self.tableWidget_8.setObjectName("tableWidget_8")
        self.tableWidget_8.setColumnCount(1)
        self.tableWidget_8.setRowCount(100)
        self.label_13 = QtWidgets.QLabel(self.tab_drill)
        self.label_13.setGeometry(QtCore.QRect(620, 30, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.tableWidget_9 = QtWidgets.QTableWidget(self.tab_drill)
        self.tableWidget_9.setGeometry(QtCore.QRect(440, 60, 531, 561))
        self.tableWidget_9.setObjectName("tableWidget_9")
        self.tableWidget_9.setColumnCount(100)
        self.tableWidget_9.setRowCount(100)
        self.tabWidget.addTab(self.tab_drill, "")

        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.tabWidget.addTab(self.tab_2, "")
        self.verticalLayout_4.addWidget(self.tabWidget)





        self.lineEdit_3 = QtWidgets.QLineEdit(self.tab)
        self.lineEdit_3.setGeometry(QtCore.QRect(210, 590, 91, 22))
        self.lineEdit_3.setObjectName("lineEdit_3")

        self.label_6 = QtWidgets.QLabel(self.tab_3)
        self.label_6.setGeometry(QtCore.QRect(30, 20, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.tableWidget_3 = QtWidgets.QTableWidget(self.tab_3)
        self.tableWidget_3.setGeometry(QtCore.QRect(20, 50, 425, 461))
        self.tableWidget_3.setObjectName("tableWidget_3")
        self.tableWidget_3.setColumnCount(4)
        # self.tableWidget_3.setRowCount(100)
        self.tableWidget_3.setHorizontalHeaderLabels(['Дата', 'От кого', 'Система координат', 'Актуальность'])

        self.label_7 = QtWidgets.QLabel(self.tab_3)
        self.label_7.setGeometry(QtCore.QRect(610, 20, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.tableWidget_4 = QtWidgets.QTableWidget(self.tab_3)
        self.tableWidget_4.setGeometry(QtCore.QRect(510, 50, 800, 461))
        self.tableWidget_4.setObjectName("tableWidget_4")
        self.tableWidget_4.setColumnCount(5)
        self.tableWidget_4.setRowCount(100)
        self.tableWidget_4.setHorizontalHeaderLabels(['Номер', 'X', 'Y', 'Пласт', 'Актуальность'])

        self.pushButton_5 = QtWidgets.QPushButton(self.tab_3)
        self.pushButton_5.setGeometry(QtCore.QRect(20, 520, 93, 28))
        self.pushButton_5.setObjectName("pushButton_5")


        self.label_2 = QtWidgets.QLabel(self.tab)
        self.label_2.setGeometry(QtCore.QRect(400, 10, 81, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.comboBox_2 = QtWidgets.QComboBox(self.tab)
        self.comboBox_2.setGeometry(QtCore.QRect(120, 50, 671, 31))
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.setEditable(True)
        self.tableView_3 = QtWidgets.QTableView(self.tab)
        self.tableView_3.setGeometry(QtCore.QRect(10, 120, 291, 461))
        self.tableView_3.setObjectName("tableView_3")
        self.pushButton_4 = QtWidgets.QPushButton(self.tab)
        self.pushButton_4.setGeometry(QtCore.QRect(10, 590, 191, 28))
        self.pushButton_4.setObjectName("pushButton_4")
        # self.pushButton_5 = QtWidgets.QPushButton(self.tab)
        # self.pushButton_5.setGeometry(QtCore.QRect(10, 630, 191, 28))
        # self.pushButton_5.setObjectName("pushButton_5")


        self.pushButton_6 = QtWidgets.QPushButton(self.tab)
        self.pushButton_6.setGeometry(QtCore.QRect(330, 590, 191, 28))
        self.pushButton_6.setObjectName("pushButton_6")

        self.pushButton_7 = QtWidgets.QPushButton(self.tab)
        self.pushButton_7.setGeometry(QtCore.QRect(10, 630, 191, 28))
        self.pushButton_7.setObjectName("pushButton_7")

        self.label_5 = QtWidgets.QLabel(self.tab)
        self.label_5.setGeometry(QtCore.QRect(400, 90, 121, 16))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")

        self.tableWidget = QtWidgets.QTableWidget(self.tab)
        self.tableWidget.setGeometry(QtCore.QRect(320, 120, 671, 461))
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(2)
        self.tableWidget.setRowCount(100)
        self.tableWidget.setHorizontalHeaderLabels(['Наименование', 'Схема'])

        self.tableWidget_2 = QtWidgets.QTableWidget(self.tab)
        self.tableWidget_2.setGeometry(QtCore.QRect(10, 120, 291, 461))
        self.tableWidget_2.setObjectName("tableWidget_2")
        self.tableWidget_2.setColumnCount(5)
        self.tableWidget_2.setRowCount(100)
        self.tableWidget_2.setHorizontalHeaderLabels(['Позиции', 'Движки', 'X', 'Y', 'id'])

        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 869, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.action = QtWidgets.QAction(MainWindow)
        self.action.setObjectName("action")
        self.action_2 = QtWidgets.QAction(MainWindow)
        self.action_2.setObjectName("action_2")
        self.action_3 = QtWidgets.QAction(MainWindow)
        self.action_3.setObjectName("action_3")

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Автоматизированное рабочее место маркшейдера"))
        self.label.setText(_translate("MainWindow", "Лицензионный участок"))
        self.pushButton.setText(_translate("MainWindow", "Обновить список лицензионных участков"))
        self.pushButton_2.setText(_translate("MainWindow", "Добавить/Изменить"))
        self.checkBox.setText(_translate("MainWindow", "Актуальность"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Wellpads), _translate("MainWindow", "Площадки"))
        self.pushButton_3.setText(_translate("MainWindow", "Изменить"))
        self.checkBox_2.setText(_translate("MainWindow", "Актуальность"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Rigs), _translate("MainWindow", "Станки"))
        self.label_2.setText(_translate("MainWindow", "Схема"))
        self.pushButton_4.setText(_translate("MainWindow", "Добавить"))
        # self.pushButton_5.setText(_translate("MainWindow", "Изменить"))
        self.pushButton_6.setText(_translate("MainWindow", "Сохранить"))

        self.pushButton_7.setText(_translate("MainWindow", "Сохранить"))


        self.label_5.setText(_translate("MainWindow", "Шаблоны схем"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "Позиции"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3), _translate("MainWindow", "Цели"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Устья"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_drill), _translate("MainWindow", "Расчет бурения"))


        self.label_6.setText(_translate("MainWindow", "Документ"))
        self.label_7.setText(_translate("MainWindow", "Цели"))
        self.pushButton_5.setText(_translate("MainWindow", "Добавить"))

        self.action.setText(_translate("MainWindow", "Добавить"))
        self.action_2.setText(_translate("MainWindow", "Добавить"))
        self.action_3.setText(_translate("MainWindow", "Добавить"))

        self.label_11.setText(_translate("MainWindow", "Позиции"))
        self.label_12.setText(_translate("MainWindow", "Цели"))
        self.label_13.setText(_translate("MainWindow", "Расчет бурения"))


class Ui_Dialog(QDialog):
    def __init__(self):
        QDialog.__init__(self)
        # Required by Qt4 to initialize the UI
        self.setupUi(self)

    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(540, 384)

        self.label_6 = QtWidgets.QLabel(Dialog)
        self.label_6.setGeometry(QtCore.QRect(30, 40, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.label_7 = QtWidgets.QLabel(Dialog)
        self.label_7.setGeometry(QtCore.QRect(30, 80, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(Dialog)
        self.label_8.setGeometry(QtCore.QRect(30, 130, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.label_9 = QtWidgets.QLabel(Dialog)
        self.label_9.setGeometry(QtCore.QRect(30, 280, 161, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.dateEdit = QtWidgets.QDateEdit(Dialog)
        self.dateEdit.setGeometry(QtCore.QRect(200, 40, 201, 22))
        self.dateEdit.setObjectName("dateEdit")

        self.dateEdit.setDate(date.today())

        self.checkBox = QtWidgets.QCheckBox(Dialog)
        self.checkBox.setGeometry(QtCore.QRect(200, 130, 81, 20))
        self.checkBox.setText("")
        self.checkBox.setObjectName("checkBox")
        self.checkBox.setChecked(True)

        self.comboBox = QtWidgets.QComboBox(Dialog)
        self.comboBox.setGeometry(QtCore.QRect(200, 80, 201, 22))
        self.comboBox.setObjectName("comboBox")
        self.comboBox.setEditable(True)

        self.comboBox_2 = QtWidgets.QComboBox(Dialog)
        self.comboBox_2.setGeometry(QtCore.QRect(200, 280, 201, 22))
        self.comboBox_2.setObjectName("comboBox_2")
        self.comboBox_2.setEditable(True)

        self.mQgsFileWidget = QgsFileWidget(Dialog)
        self.mQgsFileWidget.setGeometry(QtCore.QRect(200, 170, 171, 27))
        self.mQgsFileWidget.setObjectName("mQgsFileWidget")
        self.mQgsFileWidget_2 = QgsFileWidget(Dialog)
        self.mQgsFileWidget_2.setGeometry(QtCore.QRect(200, 220, 171, 27))
        self.mQgsFileWidget_2.setObjectName("mQgsFileWidget_2")
        self.label_10 = QtWidgets.QLabel(Dialog)
        self.label_10.setGeometry(QtCore.QRect(30, 170, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.label_11 = QtWidgets.QLabel(Dialog)
        self.label_11.setGeometry(QtCore.QRect(30, 220, 121, 21))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")

        self.pushButton = QtWidgets.QPushButton(Dialog)
        self.pushButton.setGeometry(QtCore.QRect(120, 340, 111, 28))
        self.pushButton.setObjectName("pushButton")
        self.pushButton_2 = QtWidgets.QPushButton(Dialog)
        self.pushButton_2.setGeometry(QtCore.QRect(280, 340, 111, 28))
        self.pushButton_2.setObjectName("pushButton_2")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        # Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label_6.setText(_translate("Dialog", "Дата"))
        self.label_7.setText(_translate("Dialog", "От кого"))
        self.label_8.setText(_translate("Dialog", "Актуальность"))
        self.label_9.setText(_translate("Dialog", "Система координат"))
        self.label_10.setText(_translate("Dialog", "Документ"))
        self.label_11.setText(_translate("Dialog", "Цели"))
        self.pushButton.setText(_translate("Dialog", "Сохранить"))
        self.pushButton_2.setText(_translate("Dialog", "Отмена"))

