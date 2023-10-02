import cmath

from .point import Point


class Position:
    def __init__(self, point_start: Point, dist: list, nds: float):
        self.point_start = point_start
        self.dist = dist
        self.nds = nds
        self.nds_rad = float(self.nds) * cmath.pi / 180.0

        self.positions = []

    def direction(self):
        if self.nds == 0.00000:
            return cmath.pi/2
        elif self.nds == 90:
            return 0.0
        elif self.nds == 180:
            return -cmath.pi/2
        elif self.nds == 270:
            return cmath.pi
        elif 0 < self.nds < 90:
            return cmath.pi/2 - self.nds_rad
        elif 90 < self.nds < 180:
            return -(self.nds_rad - cmath.pi/2)
        elif 180 < self.nds < 270:
            return -(self.nds_rad - cmath.pi/2)
        else:
            return 2*cmath.pi - self.nds_rad + cmath.pi/2

    def point_on_direction(self):

        posits = [self.point_start.get_point()]

        for d in self.dist:
            temp_c = cmath.rect(int(d), self.direction())
            temp_p = Point(temp_c.real, temp_c.imag)
            new_X = self.point_start.x + temp_p.x
            new_Y = self.point_start.y + temp_p.y

            new = Point(new_X, new_Y)

            self.point_start = new
            posits.append(new.get_point())
        return posits





