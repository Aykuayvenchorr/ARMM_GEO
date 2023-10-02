"""Геодезические задачи"""
import cmath

from ..modules.point import Point
# from modules.point import Point


def dir_ang_to_polar(angle: float) -> float:
    angle_rad = float(angle) * cmath.pi / 180.0
    if angle == 0.00000:
        angle = cmath.pi / 2
    elif angle == 90:
        angle = 0.0
    elif angle == 180:
        angle = -cmath.pi / 2
    elif angle == 270:
        angle = cmath.pi
    elif 0 < angle < 90:
        angle = cmath.pi / 2 - angle_rad
    elif 90 < angle < 180:
        angle = -(angle_rad - cmath.pi / 2)
    elif 180 < angle < 270:
        angle = -(angle_rad - cmath.pi / 2)
    else:
        angle = 2 * cmath.pi - angle_rad + cmath.pi / 2
    return angle


def pgz(point: Point, angle: float, dist: float) -> Point:
    angle_ = dir_ang_to_polar(angle)

    temp_c = cmath.rect(float(dist), angle_)
    temp_p = Point(temp_c.real, temp_c.imag)
    new_X = point.x + temp_p.x
    new_Y = point.y + temp_p.y

    new = Point(new_X, new_Y)

    return new


def n_point_nds(p_start: Point, nds: float, shifts: list[float]) -> list[tuple]:
    points = [p_start]
    for shift in shifts:
        p_new = pgz(p_start, nds, shift)
        points.append(p_new)
        p_start = p_new
    rig_positions = []
    for p in points:
        rig_positions.append(p.get_point())
    return rig_positions


def ogz(p1: Point, p2: Point) -> tuple[float, float]:
    # Преобразование координат точек в комплексные числа
    z1 = complex(p1.x, p1.y)
    z2 = complex(p2.x, p2.y)

    # Разность комплексных чисел для определения вектора между точками
    dz = z2 - z1

    # Расчет азимута (курса) в градусах
    azimuth = cmath.phase(dz)  # phase() возвращает аргумент комплексного числа в радианах
    azimuth_degrees = (azimuth * 180) / cmath.pi  # Преобразование в градусы

    # Расчет длины линии между точками
    distance = abs(dz)

    return azimuth_degrees, distance


# p1 = Point(0, 0)
# p2 = Point(1, 1)
# azimuth, distance = ogz(p1, p2)
# print(f"Азимут (курс): {azimuth} градусов")
# print(f"Длина линии: {distance} единиц")
