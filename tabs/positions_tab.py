import os
import uuid

from PyQt5.QtWidgets import QTableWidgetItem

plugin_dir = os.path.dirname(__file__)

file_path = os.path.join(
    plugin_dir,
    '../',
    'files',
    'scheme')


def save_table_data(tableWidget, comboBox_2):
    """Сохранить схемы в файл"""
    rows = tableWidget.rowCount()
    columns = tableWidget.columnCount()
    data = []

    for row in range(rows):
        row_data = []
        for col in range(columns):
            item = tableWidget.item(row, col)
            if item is not None:
                row_data.append(item.text())
            else:
                continue
        data.append(row_data)
    # открытие файла на дозапись, чтобы не удалялось предыдущее
    with open(file_path, 'w') as f:
        for row_data in data:
            if row_data:
                f.write('#'.join(row_data) + '\n')

    fill_scheme_cmbx(comboBox_2)


def fill_scheme_cmbx(comboBox_2):
    """Метод для заполнения комбобокса схемами"""
    scheme = []
    with open(file_path) as f:
        for line in f:
            scheme.append(line.strip().split('#')[1])

    existing_items = set(comboBox_2.itemText(i) for i in range(comboBox_2.count()))

    for schem in scheme:
        if schem not in existing_items:
            comboBox_2.addItem(schem)


# def calc_position(tableWidget_2, change_func_2):
#     """Вычисление позиций в соответствие со схемой"""
#     tableWidget_2.clear()
#     tableWidget_2.setHorizontalHeaderLabels(['Позиции', 'Движки', 'X', 'Y', 'id'])
#
#     coords_rig = self.calc_XY_positions(self.change_func_2())
#     rig_X = coords_rig[0]
#     rig_Y = coords_rig[1]
#
#     x = QTableWidgetItem(str(rig_X))
#     y = QTableWidgetItem(str(rig_Y))
#
#     self.dlg.tableWidget_2.setItem(0, 2, x)
#     self.dlg.tableWidget_2.setItem(0, 3, y)
#
#     positions = int(self.dlg.lineEdit_3.text())
#     data = []
#     i = 1
#
#     for pos in range(positions):
#         data.append(i)
#         i += 1
#     schema = ['0']
#     [schema.append(el) for el in self.dlg.comboBox_2.currentText().strip().split('_')]
#
#     if len(data) >= len(schema[1::]):
#         for row, item in enumerate(schema):
#             if row != len(schema):
#                 position = QTableWidgetItem(str(row + 1))
#                 movement = QTableWidgetItem(str(item))
#                 self.dlg.tableWidget_2.setItem(row, 0, position)
#                 self.dlg.tableWidget_2.setItem(row, 1, movement)
#
#                 self.dlg.tableWidget_2.setItem(row, 2, QTableWidgetItem(str(self.positions[row][0])))
#                 self.dlg.tableWidget_2.setItem(row, 3, QTableWidgetItem(str(self.positions[row][1])))
#                 if row != 0:
#                     self.dlg.tableWidget_2.setItem(row, 4, QTableWidgetItem(str(uuid.uuid4())))
#
#
#     else:
#         position = QTableWidgetItem('1')
#         movement = QTableWidgetItem('0')
#         self.dlg.tableWidget_2.setItem(0, 0, position)
#         self.dlg.tableWidget_2.setItem(0, 1, movement)
#
#         for row, item in enumerate(data, start=1):
#             position = QTableWidgetItem(str(item + 1))
#             movement = QTableWidgetItem(str(schema[row]))
#             self.dlg.tableWidget_2.setItem(row, 0, position)
#             self.dlg.tableWidget_2.setItem(row, 1, movement)
#             self.dlg.tableWidget_2.setItem(row, 2, QTableWidgetItem(str(self.positions[row][0])))
#             self.dlg.tableWidget_2.setItem(row, 3, QTableWidgetItem(str(self.positions[row][1])))
#             if row != 0:
#                 self.dlg.tableWidget_2.setItem(row, 4, QTableWidgetItem(str(uuid.uuid4())))
#
#     self.dlg.tableWidget_9.clearContents()
#     self.col = 0
