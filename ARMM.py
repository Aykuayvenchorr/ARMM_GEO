# -*- coding: utf-8 -*-
"""
/***************************************************************************
 ARMM
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
import uuid

from PyQt5.QtCore import QStringListModel, Qt, QAbstractTableModel, QSortFilterProxyModel
from PyQt5.QtGui import QColor, QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QComboBox, QAbstractItemView, QStyledItemDelegate, QTableView, QTableWidgetItem, QWidget, \
    QVBoxLayout, QLabel
from PyQt5.uic.properties import QtGui
from qgis.PyQt.QtCore import QSettings, QTranslator, QCoreApplication
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction
from qgis.core import (
    QgsVectorLayer,
    QgsDataSourceUri,
    QgsFeature,
    QgsProject
)

# Initialize Qt resources from file resources.py
from .modules.position import Position
from .modules.point import Point
from .modules.rig import Rig
from .modules.wellpad import Wellpad
from .modules.lic_area import LicArea
from .resources import *
# Import the code for the dialog
from .ARMM_dialog import Ui_MainWindow, Ui_Dialog
import os.path
import psycopg2
from .settings import env


class MyTableWellpadModel(QAbstractTableModel):
    def __init__(self, data, header):
        super().__init__()

        self._data = data
        self.header = header
        self.new_data = None

    def rowCount(self, parent):
        return 1

    def columnCount(self, parent):
        return len(self._data)

    def data(self, index, role):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            return str(self._data[index.column()])

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            if index.column() == 1:
                self._data[index.column()] = value

                self.dataChanged.emit(index, index)
                self.new_data = self._data
            return True
        return False

    def get_object(self):
        return self.new_data

    def flags(self, index):
        if index.column() == 1:
            return Qt.ItemIsEditable | Qt.ItemIsEnabled

        return Qt.ItemIsEnabled

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.header[section]


class MyTableRigModel(QAbstractTableModel):
    def __init__(self, data, header):
        super().__init__()

        self._data = data
        self.header = header
        self.new_data = None

    def rowCount(self, parent):
        return 1

    def columnCount(self, parent):
        return len(self._data)

    def data(self, index, role):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            return str(self._data[index.column()])

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            # if index.column() == 1:
            self._data[index.column()] = value
            return True
        return False

    def get_object(self):
        return self.new_data

    def flags(self, index):
        if index.column() in [1, 2, 4, 5, 6]:
            return Qt.ItemIsEditable | Qt.ItemIsEnabled

        return Qt.ItemIsEnabled

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.header[section]


class MyTableSchemeModel(QAbstractTableModel):
    def __init__(self, data, header):
        super().__init__()

        self._data = data
        self.header = header
        self.new_data = None

    def rowCount(self, parent):
        return len(self._data)

    def columnCount(self, parent):
        return 2

    def data(self, index, role):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            return str(self._data[index.row()])

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            # if index.column() == 1:
            self._data[index.row()] = value
            return True
        return False

    def get_object(self):
        return self.new_data

    def flags(self, index):
        if index.column() in [0, 1]:
            return Qt.ItemIsEditable | Qt.ItemIsEnabled

        return Qt.ItemIsEnabled

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.header[section]


class MyTablePositionModel(QAbstractTableModel):
    def __init__(self, data, header):
        super().__init__()

        self._data = data
        self.header = header
        self.new_data = None

    def rowCount(self, parent):
        return len(self._data)

    def columnCount(self, parent):
        return 2

    def data(self, index, role):
        if role == Qt.DisplayRole or role == Qt.EditRole:
            return str(self._data[index.column()])

    def setData(self, index, value, role):
        if role == Qt.EditRole:
            # if index.column() == 1:
            self._data[index.column()] = value
            return True
        return False

    def get_object(self):
        return self.new_data

    def flags(self, index):
        if index.column() in [0, 1]:
            return Qt.ItemIsEditable | Qt.ItemIsEnabled

        return Qt.ItemIsEnabled

    def headerData(self, section, orientation, role):
        if role == Qt.DisplayRole:
            if orientation == Qt.Horizontal:
                return self.header[section]


class ARMM:
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.positions = None
        self.rig = None
        self.dlg_2 = None
        self.cur = None
        self.conn = None
        self.text_rig = None
        self.new_text = None
        self.text_wp = None
        self.change_wp = None
        self.list_for_cmbEdit = None
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            'ARMM_{}.qm'.format(locale))

        self.file_path = os.path.join(
            self.plugin_dir,
            'files',
            'scheme')
        self.file_path_from = os.path.join(
            self.plugin_dir,
            'files',
            'from_who')
        self.file_path_crs = os.path.join(
            self.plugin_dir,
            'files',
            'crs')

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            QCoreApplication.installTranslator(self.translator)

        # Declare instance attributes
        self.actions = []
        self.menu = self.tr(u'&ARMM')
        self.active_wellpad = False

        # Check if plugin was started the first time in current QGIS session
        # Must be set in initGui() to survive plugin reloads
        self.first_start = None
        self.lic_areas = None
        self.wellpads = None
        self.rigs = None
        self.list_wp = None
        self.list_rig = None

    # noinspection PyMethodMayBeStatic

    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('ARMM', message)

    def add_action(
            self,
            icon_path,
            text,
            callback,
            enabled_flag=True,
            add_to_menu=True,
            add_to_toolbar=True,
            status_tip=None,
            whats_this=None,
            parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            # Adds plugin icon to Plugins toolbar
            self.iface.addToolBarIcon(action)

        if add_to_menu:
            self.iface.addPluginToMenu(
                self.menu,
                action)

        self.actions.append(action)

        return action

    def initGui(self):
        """Create the menu entries and toolbar icons inside the QGIS GUI."""

        icon_path = ':/plugins/ARMM/icon.png'
        self.add_action(
            icon_path,
            text=self.tr(u'ARMM'),
            callback=self.run,
            parent=self.iface.mainWindow())

        # will be set False in run()
        self.first_start = True

    def unload(self):
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginMenu(
                self.tr(u'&ARMM'),
                action)
            self.iface.removeToolBarIcon(action)

    def run(self):
        """Run method that performs all the real work"""

        # Create the dialog with elements (after translation) and keep reference
        # Only create GUI ONCE in callback, so that it will only load when the plugin is started
        if self.first_start:
            self.first_start = False
            self.dlg = Ui_MainWindow()
            self.dlg_2 = Ui_Dialog()

        self.dlg.show()

        # Установка параметров подключения к базе данных
        self.conn = psycopg2.connect(
            host="localhost",
            port="5432",
            database="arm_db",
            user="postgres",
            password="postgres"
        )
        # Создание курсора для выполнения SQL-запросов
        self.cur = self.conn.cursor()

        self.dlg.cmbLic.clear()
        self.fill_lic_areas()

        self.model_1 = QStandardItemModel()
        self.model_2 = QStandardItemModel()

        self.dlg.listView_wp.setModel(self.model_1)
        self.dlg.listView_rig.setModel(self.model_2)

        self.dlg.cmbLic.currentTextChanged.connect(self.choose_lic_area)

        self.dlg.pushButton.clicked.connect(self.fill_lic_areas)
        self.dlg.pushButton_6.clicked.connect(self.save_table_data)
        self.dlg.pushButton_4.clicked.connect(self.calc_position)

        self.dlg.listView_wp.doubleClicked.connect(self.change_func)
        self.dlg.listView_wp.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.dlg.listView_wp.doubleClicked.connect(self.get_editable_wellpad)
        self.dlg.listView_wp.doubleClicked.connect(self.fill_wellpad_table)

        self.dlg.listView_rig.doubleClicked.connect(self.get_editable_rig)
        self.dlg.listView_rig.doubleClicked.connect(self.fill_rig_table)
        self.dlg.listView_rig.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.dlg.listView_rig.doubleClicked.connect(self.change_func_2)
        self.dlg.listView_rig.doubleClicked.connect(self.zoom_to_rig)

        self.dlg.Rigs.setEnabled(False)

        self.dlg.lineEdit.textEdited.connect(self.get_new_name)
        self.dlg.lineEdit_2.textEdited.connect(self.get_new_name)

        self.dlg.pushButton_2.clicked.connect(self.edit_wellpads)
        self.dlg.pushButton_3.clicked.connect(self.edit_rigs)
        self.dlg.pushButton_5.clicked.connect(self.open_new_window)
        self.dlg.pushButton_7.clicked.connect(self.save_positions_to_BD)

        self.fill_from_cmbx()
        self.fill_crs_cmbx()
        self.dlg_2.pushButton_2.clicked.connect(self.close_dialog)
        self.dlg_2.pushButton.clicked.connect(self.save_docs)
        self.dlg_2.pushButton.clicked.connect(self.save_dialog)

        # self.dlg_2.comboBox.setEditable(True)

    def choose_lic_area(self):
        """Метод, который заполняет площадки для соответствующего ЛУ"""
        # считываем текущее значение ЛУ
        self.dlg.Rigs.setEnabled(False)
        # устанавливаем логическую переменную, которая отвечает за смену площадки (нужна, чтобы потом в станки
        # не записывались все станки, а только для одной площадки, и только один раз (один список))
        self.change_wp = True
        self.model_1.clear()
        self.model_2.clear()
        self.model_2.removeRows(0, self.model_2.columnCount())
        self.load_wellpads()

        cur_lic = self.dlg.cmbLic.currentText()
        for key, value in self.wellpads.get_dict_wellpads().items():
            # проверяем, что luid у площадки и id у ЛУ совпадают
            try:
                if value[1] == self.lic_areas.get_dict_lic_area()[cur_lic]:
                    # заносим в список все площадки, принадлежащие одному ЛУ
                    self.list_wp = self.wellpads.create_list_wellpads_for_one_lic_area(value[1])
                    if self.list_wp:
                        self.dlg.Wellpads.setEnabled(True)
                    # очищаем список с площадками, чтобы для другого ЛУ не отображались прошлые площадки
                    # self.list_wp.clear()
                    # прерываем цикл, чтобы площадки не множились
                    break
            except:
                # возникает ошибка из-за поля '-------'
                KeyError
        self.append_wellpad(self.list_wp)
        for name in self.list_wp:
            relevant = self.wellpads.get_dict_wellpads()[name][2]
            wp_id = self.wellpads.get_dict_wellpads()[name][0]
            self.unactual_rigs(relevant, wp_id)
        if self.list_wp:
            self.list_wp.clear()

        return self.list_wp

    def choose_wellpad(self, wp):
        """Метод, который заполняет станки для соответствующей плошадки"""
        self.dlg.Rigs.setEnabled(True)
        self.model_2.clear()
        self.load_rigs()
        # self.model_2.removeRows(0, self.model_2.rowCount())
        # Для того чтобы нельзя было поменять актуальность у станков, у которых площадка неактуальна
        if not self.wellpads.get_dict_wellpads()[wp][2]:
            self.dlg.checkBox_2.setEnabled(False)
        else:
            self.dlg.checkBox_2.setEnabled(True)

        try:
            for key, value in self.rigs.get_dict_rigs().items():
                if value[2] == self.wellpads.get_dict_wellpads()[wp][0]:
                    self.list_rig = self.rigs.create_list_rigs_for_one_wellpad(value[2])

                    self.change_wp = False
                    if self.list_rig:
                        self.dlg.Rigs.setEnabled(True)

                    # очищаем список с площадками, чтобы для другого ЛУ не отображались прошлые площадки
                    # self.list_rig.clear()
                    # прерываем цикл, чтобы площадки не множились
                    break
        except:
            KeyError

        self.append_rig(self.list_rig)
        if self.list_rig:
            self.list_rig.clear()
        # self.rigs.null_rigs()
        self.fill_scheme_table()
        self.fill_scheme_cmbx()
        # self.calc_XY_positions()

        # self.fill_position_table()

        return self.list_rig

    def change_func(self):
        rez = self.dlg.listView_wp.currentIndex().data()
        # self.model_2.clear()
        # self.model_2.removeRows(0, self.model_2.columnCount())
        self.choose_wellpad(rez)

        # self.model_2.removeRows(0, self.model_2.rowCount())

        return rez

    def change_func_2(self):
        rez = self.dlg.listView_rig.currentIndex().data()
        # self.calc_XY_positions(rez)
        return rez

    def update_list_lic_areas(self):
        self.dlg.cmbLic.clear()

        conn = psycopg2.connect(
            host="localhost",
            port="5432",
            database="arm_db",
            user="postgres",
            password="postgres"
        )
        cur = conn.cursor()

        # Выполнение SQL-запроса для получения всех объектов из таблицы
        cur.execute("SELECT * FROM lic_area")

        # Получение всех строк lic_area (объектов) из результата запроса
        rows = cur.fetchall()
        self.lic_areas = LicArea(rows)

        for name in self.lic_areas.get_dict_lic_area():
            self.dlg.cmbLic.addItem(name)

        cur.close()
        conn.close()

    def get_editable_wellpad(self):
        """Выбирает редактируемую площадку и устанавливает ее в lineEdit"""
        self.text_wp = self.dlg.listView_wp.currentIndex().data()
        self.dlg.lineEdit.setText(self.text_wp)
        return self.text_wp

    def get_editable_rig(self):
        """Выбирает редактируемый станок и устанавливает его в lineEdit_2"""
        self.text_rig = self.dlg.listView_rig.currentIndex().data()
        self.dlg.lineEdit_2.setText(self.text_rig)
        return self.text_rig

    def get_new_name(self, text):
        """Считывает новое имя из lineEdit и lineEdit_2"""
        self.new_text = text

    def edit_wellpads(self):
        """Метод, который позволяет добавлять и редактировать существующие площадки"""
        # conn = psycopg2.connect(
        #     host="localhost",
        #     port="5432",
        #     database="arm_db",
        #     user="postgres",
        #     password="postgres"
        # )
        # cur = conn.cursor()
        cur_lic = self.dlg.cmbLic.currentText()

        model = self.dlg.tableView_2.model()
        relevant = False

        if self.text_wp == '-------':
            self.text_wp = None
            luid = self.lic_areas.get_dict_lic_area()[cur_lic]
            if self.dlg.checkBox.isChecked():
                relevant = True
            wp = (self.new_text, luid, relevant)
            self.cur.execute("INSERT INTO wellpad (name, luid, rel) VALUES (%s, %s, %s)", wp)
            self.conn.commit()
            print("Данные добавлены")
            # self.cur.close()
            # self.conn.close()
            self.dlg.lineEdit.clear()
            # self.run()
            self.load_wellpads()
            self.dlg.cmbLic.setItemText(0, cur_lic)

        elif self.text_wp in self.wellpads.get_dict_wellpads():
            data = model._data
            new_name = self.new_text
            object_id = data[0]
            new_name_ = data[1]
            if self.dlg.checkBox.isChecked():
                relevant = True
                self.unactual_rigs(relevant, object_id)
                self.dlg.checkBox_2.setEnabled(True)
            else:
                relevant = False
                self.unactual_rigs(relevant, object_id)
                # Для того чтобы нельзя было поменять актуальность у станков, у которых площадка неактуальна
                self.dlg.checkBox_2.setEnabled(False)
            if new_name:
                sql_update_query = """Update wellpad set name = %s, rel = %s where name = %s"""
                self.cur.execute(sql_update_query, (new_name, relevant, self.text_wp))
            if new_name_ != self.text_wp:
                sql_query = """
                        UPDATE wellpad
                        SET name = %s
                        WHERE id = %s
                        """
                self.cur.execute(sql_query, (new_name_, object_id))

                self.conn.commit()
            else:
                sql_update_query = """Update wellpad set name = %s, rel = %s where name = %s"""
                self.cur.execute(sql_update_query, (self.text_wp, relevant, self.text_wp))
            self.conn.commit()
            print("Данные обновлены")
            # cur.close()
            # conn.close()
            self.dlg.lineEdit.clear()
            # перезагружаем окно, оставляя ЛУ выбранным, чтобы не выходить и заходить заново в модуль
        self.load_wellpads()
        self.choose_lic_area()
        # self.run()
        self.dlg.cmbLic.setItemText(0, cur_lic)

    def edit_rigs(self):
        """Метод, который позволяет редактировать существующие станки"""
        # conn = psycopg2.connect(
        #     host="localhost",
        #     port="5432",
        #     database="arm_db",
        #     user="postgres",
        #     password="postgres"
        # )
        # cur = conn.cursor()
        cur_lic = self.dlg.cmbLic.currentText()

        model = self.dlg.tableView.model()
        # data = model._data

        if self.text_rig in self.rigs.get_dict_rigs():
            data = model._data

            new_name = self.new_text
            relevant = False
            object_id = data[0]
            new_name_ = data[1]
            type_ = data[2]
            nds_ = data[4]
            nwells_ = data[5]
            radius_ = data[6]

            if self.dlg.checkBox_2.isChecked():
                relevant = True
            if new_name:
                sql_update_query = """Update rig set name = %s, type =  %s, nds = %s, nwells = %s, radius = %s, 
                rel = %s where name = %s"""
                self.cur.execute(sql_update_query, (new_name, type_, nds_, nwells_, radius_, relevant, self.text_rig))
            if new_name_ != self.text_rig:
                sql_query = """Update rig set name = %s, type =  %s, nds = %s, nwells = %s, radius = %s, 
                rel = %s where id = %s"""
                self.cur.execute(sql_query, (new_name_, type_, nds_, nwells_, radius_, relevant, object_id))
            else:
                sql_update_query = """Update rig set name = %s, type =  %s, nds = %s, nwells = %s, radius = %s, 
                rel = %s where name = %s"""
                self.cur.execute(sql_update_query,
                                 (self.text_rig, type_, nds_, nwells_, radius_, relevant, self.text_rig))

            self.dlg.listView_rig.setModel(self.sort_listview(self.model_2))
            self.conn.commit()
            print("Данные обновлены")
            # self.cur.close()
            # self.conn.close()
            self.dlg.lineEdit_2.clear()
            # перезагружаем окно, оставляя ЛУ выбранным, чтобы не выходить и заходить заново в модуль
            self.change_func()
            # self.run()
            self.dlg.cmbLic.setItemText(0, cur_lic)
            # как сделать так, чтобы не надо было заново заходить в площадки
            # self.dlg.listView_wp.setItem(self.text_wp)

    def fill_wellpad_table(self):
        """Метод для детального отображения площадки в tableview"""
        rez = self.text_wp
        header = ['№', 'Название', 'ID ЛУ', 'Актуальность']
        try:
            data = [self.wellpads.get_dict_wellpads()[rez][0], rez, self.wellpads.get_dict_wellpads()[rez][1],
                    self.wellpads.get_dict_wellpads()[rez][2]]

            model = MyTableWellpadModel(data, header)

            self.dlg.tableView_2.setModel(model)
        except:
            KeyError

    def fill_rig_table(self):
        """Метод для детального отображения станка в tableview"""
        rez = self.text_rig
        header = ['№', 'Название', 'Тип', 'ID площадки', 'НДС', 'Радиус', 'Актуальность']
        try:
            data = [self.rigs.get_dict_rigs()[rez][0], rez, self.rigs.get_dict_rigs()[rez][1],
                    self.rigs.get_dict_rigs()[rez][2], self.rigs.get_dict_rigs()[rez][3],
                    self.rigs.get_dict_rigs()[rez][5],
                    self.rigs.get_dict_rigs()[rez][6]]

            model = MyTableRigModel(data, header)

            self.dlg.tableView.setModel(model)
        except:
            KeyError

    # def fill_position_table(self):
    #     """Метод для детального отображения станка в tableview"""
    #     header = ['Позиции', 'Движки']
    #     try:
    #         data = [1, 2]
    #
    #         model = MyTableSchemeModel(data, header)
    #
    #         self.dlg.tableView_3.setModel(model)
    #     except:
    #         KeyError

    def fill_scheme_table(self):
        """Метод для заполнения таблицы шаблона схем"""
        # self.dlg.tableWidget.clear()
        self.dlg.tableWidget.setHorizontalHeaderLabels(['Наименование', 'Схема'])

        data = []
        with open(self.file_path) as file:
            for line in file:
                line = line.strip().split('#')
                # print(line)
                data.append((line[0], line[1]))
                # data = [
                #     ('нулевая', '0_0_0'),
                #     ('5_на_15', '5_5_5_15'),
                # ]
        for row, item in enumerate(data):
            name_item = QTableWidgetItem(item[0])
            schema_item = QTableWidgetItem(item[1])
            self.dlg.tableWidget.setItem(row, 0, name_item)
            self.dlg.tableWidget.setItem(row, 1, schema_item)

    def append_wellpad(self, list_wellpads):
        """Добавление площадок в listview и закрашивание их в соответствии с актуальностью"""
        try:
            for el in list_wellpads:
                item = QStandardItem(el)
                if not self.wellpads.get_dict_wellpads()[el][2]:
                    item.setForeground(QColor(Qt.gray))
                    self.model_1.appendRow(item)
                else:
                    self.model_1.appendRow(item)
        except:
            KeyError
        item = QStandardItem('-------')
        self.model_1.appendRow(item)

        self.dlg.listView_wp.setModel(self.sort_listview(self.model_1))

    def append_rig(self, list_rigs):
        """Добавление станков в listview и закрашивание их в соответствии с актуальностью"""
        try:
            for el in list_rigs:
                item = QStandardItem(el)
                if not self.rigs.get_dict_rigs()[el][6]:
                    item.setForeground(QColor(Qt.gray))
                    self.model_2.appendRow(item)
                else:
                    self.model_2.appendRow(item)
        except:
            KeyError
        item = QStandardItem('-------')
        self.model_2.appendRow(item)

        self.dlg.listView_rig.setModel(self.sort_listview(self.model_2))

    def sort_listview(self, model):
        """Метод для сортировки названий по алфавиту"""
        proxy_model = QSortFilterProxyModel()
        proxy_model.setSourceModel(model)
        proxy_model.setSortCaseSensitivity(Qt.CaseInsensitive)
        proxy_model.sort(0)
        return proxy_model

    def zoom_to_rig(self):
        self.iface.mapCanvas().refresh()
        cur_rig = self.change_func_2()
        cur_rig_id = self.rigs.get_dict_rigs()[cur_rig][0]
        # project_path = "D:\ГПН_ГЕО\ARMM\\rtmm_test.qgz"

        # Загрузка проекта QGIS
        project = QgsProject.instance()
        # project.read(project_path)

        # Получение слоя по имени "rig"
        layer_name = "rig"
        layer = QgsProject.instance().mapLayersByName(layer_name)[0]
        # Снимаем выделения со всех станков, чтобы выделялся только один
        layer.removeSelection()

        # Проверка, что слой успешно получен
        if layer is not None:
            layer.select([cur_rig_id])
            self.iface.mapCanvas().zoomToFeatureIds(layer, [cur_rig_id])
        else:
            print(f"Слой {layer_name} не найден!")

    def load_wellpads(self):
        """Метод для получения площадок из БД"""
        self.cur.execute("SELECT * FROM wellpad")
        rows = self.cur.fetchall()
        self.wellpads = Wellpad(rows)

    def load_rigs(self):
        """Метод для получения станков из БД"""
        self.cur.execute("SELECT * FROM rig")
        rows = self.cur.fetchall()
        self.rigs = Rig(rows)
        # self.cur.close()
        # self.conn.close()

    def fill_lic_areas(self):
        """Метод для получения лицензионных участков из БД"""
        # Выполнение SQL-запроса для получения всех объектов из таблицы
        self.dlg.cmbLic.clear()
        self.cur.execute("SELECT * FROM lic_area")

        # Получение всех строк lic_area (объектов) из результата запроса
        rows = self.cur.fetchall()
        self.lic_areas = LicArea(rows)

        for name in self.lic_areas.get_dict_lic_area():
            self.dlg.cmbLic.addItem(name)

    def unactual_rigs(self, relevant, wp_id):
        sql_update_query = """Update rig set rel = %s where wellpad_id = %s"""
        self.cur.execute(sql_update_query, (relevant, wp_id))
        # if not relevant:
        #     self.dlg.checkBox_2.setEnabled(False)
        # if relevant:
        #     self.dlg.checkBox_2.setEnabled(True)
        # elif relevant:
        #     self.dlg.checkBox_2.setEnabled(True)

    def save_table_data(self):
        """Сохранить схемы в файл"""
        rows = self.dlg.tableWidget.rowCount()
        columns = self.dlg.tableWidget.columnCount()
        data = []

        for row in range(rows):
            row_data = []
            for col in range(columns):
                item = self.dlg.tableWidget.item(row, col)
                if item is not None:
                    row_data.append(item.text())
                else:
                    # row_data.append('')
                    continue
            data.append(row_data)
        # открытие файла на дозапись, чтобы не удалялось предыдущее
        with open(self.file_path, 'w') as f:
            for row_data in data:
                if row_data:
                    f.write('#'.join(row_data) + '\n')

        self.fill_scheme_cmbx()

    def fill_scheme_cmbx(self):
        """Метод для заполнения комбобокса схемами"""
        scheme = []
        with open(self.file_path) as f:
            for line in f:
                scheme.append(line.strip().split('#')[1])

        existing_items = set(self.dlg.comboBox_2.itemText(i) for i in range(self.dlg.comboBox_2.count()))

        for schem in scheme:
            if schem not in existing_items:
                self.dlg.comboBox_2.addItem(schem)

    def calc_position(self):
        """Вычисление позиций в соответствие со схемой"""
        self.dlg.tableWidget_2.clear()
        self.dlg.tableWidget_2.setHorizontalHeaderLabels(['Позиции', 'Движки', 'X', 'Y', 'id'])

        coords_rig = self.calc_XY_positions(self.change_func_2())
        rig_X = coords_rig[0]
        rig_Y = coords_rig[1]

        x = QTableWidgetItem(str(rig_X))
        y = QTableWidgetItem(str(rig_Y))

        self.dlg.tableWidget_2.setItem(0, 2, x)
        self.dlg.tableWidget_2.setItem(0, 3, y)

        # print(x)
        positions = int(self.dlg.lineEdit_3.text())
        data = []
        i = 1

        for pos in range(positions):
            data.append(i)
            i += 1
        schema = ['0']
        [schema.append(el) for el in self.dlg.comboBox_2.currentText().strip().split('_')]
        # print(schema)

        if len(data) >= len(schema[1::]):
            for row, item in enumerate(schema):
                if row != len(schema):
                    position = QTableWidgetItem(str(row + 1))
                    movement = QTableWidgetItem(str(item))
                    self.dlg.tableWidget_2.setItem(row, 0, position)
                    self.dlg.tableWidget_2.setItem(row, 1, movement)

                    self.dlg.tableWidget_2.setItem(row, 2, QTableWidgetItem(str(self.positions[row][0])))
                    self.dlg.tableWidget_2.setItem(row, 3, QTableWidgetItem(str(self.positions[row][1])))
                    if row != 0:
                        self.dlg.tableWidget_2.setItem(row, 4, QTableWidgetItem(str(uuid.uuid4())))

            # self.dlg.tableWidget_2.setItem(0, 2, x)
            # self.dlg.tableWidget_2.setItem(0, 3, y)

        else:
            position = QTableWidgetItem('1')
            movement = QTableWidgetItem('0')
            self.dlg.tableWidget_2.setItem(0, 0, position)
            self.dlg.tableWidget_2.setItem(0, 1, movement)
            # self.dlg.tableWidget_2.setItem(row, 2, QTableWidgetItem(str(self.positions[row][0])))
            # self.dlg.tableWidget_2.setItem(row, 3, QTableWidgetItem(str(self.positions[row][1])))
            for row, item in enumerate(data, start=1):
                position = QTableWidgetItem(str(item + 1))
                movement = QTableWidgetItem(str(schema[row]))
                self.dlg.tableWidget_2.setItem(row, 0, position)
                self.dlg.tableWidget_2.setItem(row, 1, movement)
                self.dlg.tableWidget_2.setItem(row, 2, QTableWidgetItem(str(self.positions[row][0])))
                self.dlg.tableWidget_2.setItem(row, 3, QTableWidgetItem(str(self.positions[row][1])))
                if row != 0:
                    self.dlg.tableWidget_2.setItem(row, 4, QTableWidgetItem(str(uuid.uuid4())))

            # self.dlg.tableWidget_2.setItem(0, 2, x)
            # self.dlg.tableWidget_2.setItem(0, 3, y)

    def open_new_window(self):
        self.dlg_2.exec_()

    def fill_from_cmbx(self):
        """Метод для заполнения комбобокса от кого в диалоге"""

        names = []
        with open(self.file_path_from) as f:
            for line in f:
                names.append(line.strip())

        existing_items = set(self.dlg_2.comboBox.itemText(i) for i in range(self.dlg_2.comboBox.count()))

        for name in names:
            # if schem not in existing_items:
            self.dlg_2.comboBox.addItem(name)

    def fill_crs_cmbx(self):
        """Метод для заполнения комбобокса систем координат в диалоге"""

        crs = []
        with open(self.file_path_crs) as f:
            for line in f:
                crs.append(line.strip())

        existing_items = set(self.dlg_2.comboBox_2.itemText(i) for i in range(self.dlg_2.comboBox_2.count()))

        for cr in crs:
            # if schem not in existing_items:
            self.dlg_2.comboBox_2.addItem(cr)

    def close_dialog(self):
        self.dlg_2.close()

    def save_docs(self):
        pass

    def save_dialog(self):
        """Сохраняет информацию о документе в таблицу в поле Цели"""

        date_ = QTableWidgetItem(str(self.dlg_2.dateEdit.date().toPyDate()))
        from_who = QTableWidgetItem(str(self.dlg_2.comboBox.currentText()))
        rel = QTableWidgetItem(str(self.dlg_2.checkBox.isChecked()))
        self.dlg.tableWidget_3.insertRow(0)
        self.dlg.tableWidget_3.setItem(0, 0, date_)
        self.dlg.tableWidget_3.setItem(0, 1, from_who)
        self.dlg.tableWidget_3.setItem(0, 3, rel)

    def calc_XY_positions(self, rig):
        """Рассчитывает координаты позиций станка"""
        # rig = self.change_func_2()
        rig = str(rig)
        # print(rig)

        # self.cur.execute(f"SELECT ST_AsText(geom) FROM rig WHERE name='{rig}'")
        self.cur.execute(f"SELECT ST_X(geom), ST_Y(geom) FROM rig WHERE name='{rig}'")
        rows = self.cur.fetchall()
        # self.cur.execute(f"SELECT nds FROM rig WHERE name='{rig}'")
        # nds = float(self.cur.fetchall()[0][0])
        nds = self.rigs.get_dict_rigs()[rig][3]
        # print(nds)

        # self.rig = Rig(rows)
        x = rows[0][0]
        y = rows[0][1]
        rig_pos_0 = Point(x, y)

        positions = int(self.dlg.lineEdit_3.text())

        schema = self.dlg.comboBox_2.currentText().strip().split('_')[0:positions]

        calc = Position(rig_pos_0, schema, nds)

        # print(calc.point_on_direction())

        self.positions = calc.point_on_direction()

        # print(x)
        # print(y)
        # print(rig_pos_0.x)
        return x, y
        # print(self.rig.get_dict_rigs()[rig][-1])
        # print(x)
        # print(y)

    def save_positions_to_BD(self):
        """Сохраняет позиции станка в БД"""
        rows = self.dlg.tableWidget_2.rowCount()
        columns = self.dlg.tableWidget_2.columnCount()
        data = []

        for row in range(1, rows):
            row_data = []
            item_row = self.dlg.tableWidget_2.item(row, 0)
            if item_row is not None:

                for col in range(columns):
                    item = self.dlg.tableWidget_2.item(row, col)
                    if item is not None:
                        row_data.append(item.text())
                    else:
                        # row_data.append('')
                        continue
                data.append(row_data)
        rig_id_ = self.rigs.get_dict_rigs()[self.change_func_2()][0]
        print(data)
        positions = [1]

        for pos in data:
            self.cur.execute(
                f"INSERT INTO position (number, geom, rig_id) VALUES ('{pos[0]}', 'Point({float(pos[2])} {float(pos[3])})', '{rig_id_}')")
            self.conn.commit()
            print("Данные добавлены")

            positions.append(pos[0])

        self.fill_pos_in_calc_drill(positions)

    def fill_pos_in_calc_drill(self, positions):
        """Заполняет таблицу с позициями в поле Расчет бурения"""
        # self.dlg.tableWidget_7.insertRow(0)

        for row, item in enumerate(positions):
            position = QTableWidgetItem(str(item))
            self.dlg.tableWidget_7.setItem(row, 0, position)
