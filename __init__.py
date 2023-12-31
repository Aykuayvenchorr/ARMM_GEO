# -*- coding: utf-8 -*-
"""
/***************************************************************************
 ARMM
                                 A QGIS plugin
 ARMM
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                             -------------------
        begin                : 2023-06-29
        copyright            : (C) 2023 by GPN_GEO
        email                : nastyashevchenkomail@mail.ru
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load ARMM class from file ARMM.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .ARMM import ARMM
    return ARMM(iface)
