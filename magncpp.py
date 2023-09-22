import math

gh1 = [0.0] * MAXCOEFF
gh2 = [0.0] * MAXCOEFF
gha = [0.0] * MAXCOEFF
ghb = [0.0] * MAXCOEFF

mdfile = ""
need_to_read_model = 1
model = [""] * MAXMOD
nmodel = 0
max1 = [0] * MAXMOD
max2 = [0] * MAXMOD
max3 = [0] * MAXMOD
nmax = 0
irec_pos = [0] * MAXMOD
epoch = [0.0] * MAXMOD
yrmin = [0.0] * MAXMOD
yrmax = [0.0] * MAXMOD
minyr = 0.0
maxyr = 0.0
altmin = [0.0] * MAXMOD
altmax = [0.0] * MAXMOD
minalt = 0.0
maxalt = 0.0
igdgc = 3
warn_H = 0
warn_H_strong = 0
warn_P = 0
warn_H_val = 0.0
warn_H_strong_val = 0.0
"""
/****************************************************************************/
/*                                                                          */
/*                       Subroutine degrees_to_decimal                      */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Converts degrees,minutes, seconds to decimal degrees.                */
/*                                                                          */
/*     Input:                                                               */
/*            degrees - Integer degrees                                     */
/*            minutes - Integer minutes                                     */
/*            seconds - Integer seconds                                     */
/*                                                                          */
/*     Output:                                                              */
/*            decimal - degrees in decimal degrees                          */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 12, 1988                                                */
/*                                                                          */
/****************************************************************************/
"""


def sqrmy2(x):
    return x * x


def degrees_to_decimal(degrees, minutes, seconds):
    deg = float(degrees)
    min = float(minutes) / 60.0
    sec = float(seconds) / 3600.0
    decimal = abs(sec) + abs(min) + abs(deg)

    if deg < 0:
        decimal = -decimal
    elif deg == 0:
        if min < 0:
            decimal = -decimal
        elif min == 0:
            if sec < 0:
                decimal = -decimal

    return decimal


"""
/****************************************************************************/
/*                                                                          */
/*                           Subroutine julday                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the decimal day of year from month, day, year.              */
/*     Supplied by Daniel Bergstrom                                         */
/*                                                                          */
/* References:                                                              */
/*                                                                          */
/* 1. Nachum Dershowitz and Edward M. Reingold, Calendrical Calculations,   */
/*    Cambridge University Press, 3rd edition, ISBN 978-0-521-88540-9.      */
/*                                                                          */
/* 2. Claus Tøndering, Frequently Asked Questions about Calendars,          */
/*    Version 2.9, http://www.tondering.dk/claus/calendar.html              */
/*                                                                          */
/****************************************************************************/
"""


def julday(month, day, year):
    days = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]

    leap_year = (year % 4 == 0) and ((year % 100 != 0) or (year % 400 == 0))

    day_in_year = days[month - 1] + day + (leap_year if month > 2 else 0)

    return year + day_in_year / (365.0 + leap_year)


"""
/****************************************************************************/
/*                                                                          */
/*                           Subroutine getshc                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Reads spherical harmonic coefficients from the specified             */
/*     model into an array.                                                 */
/*                                                                          */
/*     Input:                                                               */
/*           stream     - Logical unit number                               */
/*           iflag      - Flag for SV equal to ) or not equal to 0          */
/*                        for designated read statements                    */
/*           strec      - Starting record number to read from model         */
/*           nmax_of_gh - Maximum degree and order of model                 */
/*                                                                          */
/*     Output:                                                              */
/*           gh1 or 2   - Schmidt quasi-normal internal spherical           */
/*                        harmonic coefficients                             */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Bill Flanagan                                                  */
/*           NOAA CORPS, DESDIS, NGDC, 325 Broadway, Boulder CO.  80301     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 15, 1988                                                */
/*                                                                          */
/****************************************************************************/
"""


def getshc(file, iflag, strec, nmax_of_gh, gh):
    inbuff = ""
    irat = ""
    ii = 0
    ios = 0
    line_num = 0
    g = 0.0
    hh = 0.0
    trash = 0.0

    try:
        with open(file, "rt") as stream:
            stream.seek(strec)

            for nn in range(1, nmax_of_gh + 1):
                for mm in range(nn + 1):
                    if iflag == 1:
                        inbuff = stream.readline()
                        n, m, g, hh, trash, trash, irat, line_num = map(str.split, inbuff.split())
                        n, m, g, hh, line_num = int(n), int(m), float(g), float(hh), int(line_num)
                    else:
                        inbuff = stream.readline()
                        n, m, trash, trash, g, hh, irat, line_num = map(str.split, inbuff.split())
                        n, m, g, hh, line_num = int(n), int(m), float(g), float(hh), int(line_num)

                    if nn != n or mm != m:
                        ios = -2
                        return ios

                    ii += 1
                    if gh == 1:
                        gh1[ii] = g
                    elif gh == 2:
                        gh2[ii] = g
                    else:
                        print("\nError in subroutine getshc")

                    if m != 0:
                        ii += 1
                        if gh == 1:
                            gh1[ii] = hh
                        elif gh == 2:
                            gh2[ii] = hh
                        else:
                            print("\nError in subroutine getshc")

    except FileNotFoundError:
        print("\nError on opening file", file)

    return ios


"""
/****************************************************************************/
/*                                                                          */
/*                           Subroutine extrapsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Extrapolates linearly a spherical harmonic model with a              */
/*     rate-of-change model.                                                */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of base model                                  */
/*           nmax1    - maximum degree and order of base model              */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of base model                 */
/*           nmax2    - maximum degree and order of rate-of-change model    */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of rate-of-change model       */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - Schmidt quasi-normal internal spherical             */
/*                    harmonic coefficients                                 */
/*           nmax   - maximum degree and order of resulting model           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 16, 1988                                                */
/*                                                                          */
/****************************************************************************/
"""


def extrapsh(date, dte1, nmax1, nmax2, gh):
    factor = date - dte1
    if nmax1 == nmax2:
        k = nmax1 * (nmax1 + 2)
        nmax = nmax1
    else:
        if nmax1 > nmax2:
            k = nmax2 * (nmax2 + 2)
            l = nmax1 * (nmax1 + 2)
            if gh == 3:
                for ii in range(k + 1, l + 1):
                    gha[ii] = gh1[ii]
            elif gh == 4:
                for ii in range(k + 1, l + 1):
                    ghb[ii] = gh1[ii]
            else:
                print("\nError in subroutine extrapsh")
            nmax = nmax1
        else:
            k = nmax1 * (nmax1 + 2)
            l = nmax2 * (nmax2 + 2)
            if gh == 3:
                for ii in range(k + 1, l + 1):
                    gha[ii] = factor * gh2[ii]
            elif gh == 4:
                for ii in range(k + 1, l + 1):
                    ghb[ii] = factor * gh2[ii]
            else:
                print("\nError in subroutine extrapsh")
            nmax = nmax2

    if gh == 3:
        for ii in range(1, k + 1):
            gha[ii] = gh1[ii] + factor * gh2[ii]
    elif gh == 4:
        for ii in range(1, k + 1):
            ghb[ii] = gh1[ii] + factor * gh2[ii]
    else:
        print("\nError in subroutine extrapsh")

    return nmax


"""
/****************************************************************************/
/*                                                                          */
/*                           Subroutine interpsh                            */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Interpolates linearly, in time, between two spherical harmonic       */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           date     - date of resulting model (in decimal year)           */
/*           dte1     - date of earlier model                               */
/*           nmax1    - maximum degree and order of earlier model           */
/*           gh1      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of earlier model              */
/*           dte2     - date of later model                                 */
/*           nmax2    - maximum degree and order of later model             */
/*           gh2      - Schmidt quasi-normal internal spherical             */
/*                      harmonic coefficients of internal model             */
/*                                                                          */
/*     Output:                                                              */
/*           gha or b - coefficients of resulting model                     */
/*           nmax     - maximum degree and order of resulting model         */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/
"""


def interpsh(date, dte1, nmax1, dte2, nmax2, gh):
    factor = (date - dte1) / (dte2 - dte1)
    if nmax1 == nmax2:
        k = nmax1 * (nmax1 + 2)
        nmax = nmax1
    else:
        if nmax1 > nmax2:
            k = nmax2 * (nmax2 + 2)
            l = nmax1 * (nmax1 + 2)
            if gh == 3:
                for ii in range(k + 1, l + 1):
                    gha[ii] = gh1[ii] + factor * (-gh1[ii])
            elif gh == 4:
                for ii in range(k + 1, l + 1):
                    ghb[ii] = gh1[ii] + factor * (-gh1[ii])
            else:
                print("\nError in subroutine interpsh")
            nmax = nmax1
        else:
            k = nmax1 * (nmax1 + 2)
            l = nmax2 * (nmax2 + 2)
            if gh == 3:
                for ii in range(k + 1, l + 1):
                    gha[ii] = factor * gh2[ii]
            elif gh == 4:
                for ii in range(k + 1, l + 1):
                    ghb[ii] = factor * gh2[ii]
            else:
                print("\nError in subroutine interpsh")
            nmax = nmax2

    if gh == 3:
        for ii in range(1, k + 1):
            gha[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii])
    elif gh == 4:
        for ii in range(1, k + 1):
            ghb[ii] = gh1[ii] + factor * (gh2[ii] - gh1[ii])
    else:
        print("\nError in subroutine interpsh")

    return nmax


"""
/****************************************************************************/
/*                                                                          */
/*                           Subroutine shval3                              */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Calculates field components from spherical harmonic (sh)             */
/*     models.                                                              */
/*                                                                          */
/*     Input:                                                               */
/*           igdgc     - indicates coordinate system used; set equal        */
/*                       to 1 if geodetic, 2 if geocentric                  */
/*           latitude  - north latitude, in degrees                         */
/*           longitude - east longitude, in degrees                         */
/*           elev      - WGS84 altitude above ellipsoid (igdgc=1), or       */
/*                       radial distance from earth's center (igdgc=2)      */
/*           a2,b2     - squares of semi-major and semi-minor axes of       */
/*                       the reference spheroid used for transforming       */
/*                       between geodetic and geocentric coordinates        */
/*                       or components                                      */
/*           nmax      - maximum degree and order of coefficients           */
/*           iext      - external coefficients flag (=0 if none)            */
/*           ext1,2,3  - the three 1st-degree external coefficients         */
/*                       (not used if iext = 0)                             */
/*                                                                          */
/*     Output:                                                              */
/*           x         - northward component                                */
/*           y         - eastward component                                 */
/*           z         - vertically-downward component                      */
/*                                                                          */
/*     based on subroutine 'igrf' by D. R. Barraclough and S. R. C. Malin,  */
/*     report no. 71/1, institute of geological sciences, U.K.              */
/*                                                                          */
/*     FORTRAN                                                              */
/*           Norman W. Peddie                                               */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 17, 1988                                                */
/*                                                                          */
/****************************************************************************/
"""

import math


def shval3(igdgc, flat, flon, elev, nmax, gh, iext, ext1, ext2, ext3):
    earths_radius = 6371.2
    dtr = 0.01745329
    slat = math.sin(flat * dtr)

    if abs(90.0 - flat) < 0.001:
        aa = 89.999  # 300 ft. from North pole
    elif abs(90.0 + flat) < 0.001:
        aa = -89.999  # 300 ft. from South pole
    else:
        aa = flat

    clat = math.cos(aa * dtr)
    sl = [0.0] * 14
    cl = [0.0] * 14
    p = [0.0] * 119
    q = [0.0] * 119

    x = 0.0
    y = 0.0
    z = 0.0
    xtemp = 0.0
    ytemp = 0.0
    ztemp = 0.0
    sd = 0.0
    cd = 1.0
    l = 1
    n = 0
    m = 1
    npq = (nmax * (nmax + 3)) // 2

    if igdgc == 1:
        aa = a2 * clat * clat
        bb = b2 * slat * slat
        cc = aa + bb
        dd = math.sqrt(cc)
        argument = elev * (elev + 2.0 * dd) + (a2 * aa + b2 * bb) / cc
        r = math.sqrt(argument)
        cd = (elev + dd) / r
        sd = (a2 - b2) / dd * slat * clat / r

    ratio = earths_radius / r
    aa = math.sqrt(3.0)
    p[1] = 2.0 * slat
    p[2] = 2.0 * clat
    p[3] = 4.5 * slat * slat - 1.5
    p[4] = 3.0 * aa * clat * slat
    q[1] = -clat
    q[2] = slat
    q[3] = -3.0 * clat * slat
    q[4] = aa * (slat * slat - clat * clat)

    for k in range(1, npq + 1):
        if n < m:
            m = 0
            n = n + 1
            rr = ratio ** (n + 2)
            fn = float(n)

        fm = float(m)

        if k >= 5:
            if m == n:
                aa = math.sqrt(1.0 - 0.5 / fm)
                j = k - n - 1
                p[k] = (1.0 + 1.0 / fm) * aa * clat * p[j]
                q[k] = aa * (clat * q[j] + slat / fm * p[j])
                sl[m] = sl[m - 1] * cl[1] + cl[m - 1] * sl[1]
                cl[m] = cl[m - 1] * cl[1] - sl[m - 1] * sl[1]
            else:
                argument = fn ** 2 - fm ** 2
                aa = math.sqrt(argument)
                argument = ((fn - 1.0) ** 2) - (fm ** 2)
                bb = math.sqrt(argument) / aa
                cc = (2.0 * fn - 1.0) / aa
                ii = k - n
                j = k - 2 * n + 1
                p[k] = (fn + 1.0) * (cc * slat / fn * p[ii] - bb / (fn - 1.0) * p[j])
                q[k] = cc * (slat * q[ii] - clat / fn * p[ii]) - bb * q[j]

        aa = rr * gha[l]

        if m == 0:
            x = x + aa * q[k]
            z = z - aa * p[k]
            l = l + 1
        else:
            bb = rr * gha[l + 1]
            cc = aa * cl[m] + bb * sl[m]
            x = x + cc * q[k]
            z = z - cc * p[k]

            if clat > 0:
                y = y + (aa * sl[m] - bb * cl[m]) * fm * p[k] / ((fn + 1.0) * clat)
            else:
                y = y + (aa * sl[m] - bb * cl[m]) * q[k] * slat
            l = l + 2

    if iext != 0:
        aa = ext2 * cl[1] + ext3 * sl[1]
        x = x - ext1 * clat + aa * slat
        y = y + ext2 * sl[1] - ext3 * cl[1]
        z = z + ext1 * slat + aa * clat

    aa = x
    x = x * cd + z * sd
    z = z * cd - aa * sd

    return x, y, z, xtemp, ytemp, ztemp


"""
/****************************************************************************/
/*                                                                          */
/*                           Subroutine dihf                                */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*     Computes the geomagnetic d, i, h, and f from x, y, and z.            */
/*                                                                          */
/*     Input:                                                               */
/*           x  - northward component                                       */
/*           y  - eastward component                                        */
/*           z  - vertically-downward component                             */
/*                                                                          */
/*     Output:                                                              */
/*           d  - declination                                               */
/*           i  - inclination                                               */
/*           h  - horizontal intensity                                      */
/*           f  - total intensity                                           */
/*                                                                          */
/*     FORTRAN                                                              */
/*           A. Zunde                                                       */
/*           USGS, MS 964, box 25046 Federal Center, Denver, CO.  80225     */
/*                                                                          */
/*     C                                                                    */
/*           C. H. Shaffer                                                  */
/*           Lockheed Missiles and Space Company, Sunnyvale CA              */
/*           August 22, 1988                                                */
/*                                                                          */
/****************************************************************************/
"""

import math


def dihf(gh, x, y, z, xtemp, ytemp, ztemp):
    ios = gh
    sn = 0.0001
    d = None
    f = None
    h = None
    i = None
    dtemp = None
    ftemp = None
    htemp = None
    itemp = None

    for j in range(1, 2):
        if gh == 3:
            h2 = x * x + y * y
            h = math.sqrt(h2)
            f2 = h2 + z * z
            f = math.sqrt(f2)

            if f < sn:
                d = math.nan
                i = math.nan
            else:
                i = math.atan2(z, h)

                if h < sn:
                    d = math.nan
                else:
                    hpx = h + x
                    if hpx < sn:
                        d = math.pi
                    else:
                        d = 2.0 * math.atan2(y, hpx)

        elif gh == 4:
            h2 = xtemp * xtemp + ytemp * ytemp
            htemp = math.sqrt(h2)
            f2 = h2 + ztemp * ztemp
            ftemp = math.sqrt(f2)

            if ftemp < sn:
                dtemp = math.nan
                itemp = math.nan
            else:
                itemp = math.atan2(ztemp, htemp)

                if htemp < sn:
                    dtemp = math.nan
                else:
                    hpx = htemp + xtemp
                    if hpx < sn:
                        dtemp = math.pi
                    else:
                        dtemp = 2.0 * math.atan2(ytemp, hpx)

    return ios, d, f, h, i, dtemp, ftemp, htemp, itemp


def read_model(file_name):
    try:
        with open(file_name, 'r') as stream:
            fileline = 0  # First line will be 1
            modelI = -1  # First model will be 0
            irec_pos = []
            model = []
            epoch = []
            max1 = []
            max2 = []
            max3 = []
            yrmin = []
            yrmax = []
            altmin = []
            altmax = []
            minyr = None
            maxyr = None

            for line in stream:
                fileline += 1  # On new line

                if len(line) != RECL:  # If incorrect record size
                    print(f"Corrupt record in file {file_name} on line {fileline}. Need ,read {len(line)}")
                    return -5

                if line.startswith("   "):  # If 1st 3 chars are spaces
                    modelI += 1  # New model

                    if modelI > MAXMOD:  # If too many headers
                        print(f"Too many models in file {file_name} on line {fileline}.")
                        return -6

                    irec_pos.append(stream.tell())
                    fields = line.split()
                    model.append(fields[0])
                    epoch.append(float(fields[1]))
                    max1.append(int(fields[2]))
                    max2.append(int(fields[3]))
                    max3.append(int(fields[4]))
                    yrmin.append(float(fields[5]))
                    yrmax.append(float(fields[6]))
                    altmin.append(float(fields[7]))
                    altmax.append(float(fields[8]))

                    # Compute date range for all models
                    if modelI == 0:  # If first model
                        minyr = yrmin[0]
                        maxyr = yrmax[0]
                    else:
                        if yrmin[modelI] < minyr:
                            minyr = yrmin[modelI]
                        if yrmax[modelI] > maxyr:
                            maxyr = yrmax[modelI]

            nmodel = modelI + 1
            need_to_read_model = False
            return 1

    except FileNotFoundError:
        print(f"Error: File {file_name} not found.")
        return -1


# Constants
RECL = 80  # Record length
MAXMOD = 100  # Maximum number of models


def calc_point1(mdfilein, sdate, igdgc, latitude, longitude, alt, x, y, z):
    x[0] = sdate
    y[0] = latitude
    z[0] = longitude


import ctypes


# Структура для хранения моделей
class Model:
    def __init__(self):
        self.model = ""
        self.epoch = 0.0
        self.max1 = 0
        self.max2 = 0
        self.max3 = 0
        self.yrmin = 0.0
        self.yrmax = 0.0
        self.altmin = 0.0
        self.altmax = 0.0

    # Заглушки для необходимых функций
    def read_model(stream):
        # Реализуйте чтение моделей из файла
        pass

    def calculate(sdate, modelI, igdgc, latitude, longitude, alt, x, y, z, d, f, h, i, ddot, fdot, hdot, idot, xdot,
                  ydot,
                  zdot):
        # Реализуйте расчет магнитных компонентов
        pass

    # Заглушка для strncpy
    def strncpy(dest, source, count):
        # Реализуйте копирование строки
        pass


# Определение констант и переменных
MAXMOD = 100
mdfile = ""  # Путь к файлу модели
need_to_read_model = True
models = [Model() for _ in range(MAXMOD)]
nmodel = 0
minalt = 0.0
maxalt = 0.0
igdgc = 3
warn_H = 0
warn_H_val = 99999.0
warn_H_strong = 0
warn_H_strong_val = 99999.0
warn_P = 0

# Открытие файла модели
stream = None

# Выбор модели
for modelI in range(nmodel):
    if sdate < models[modelI].yrmax:
        break

if modelI == nmodel:
    modelI -= 1

# Получение минимальной и максимальной высоты для выбранной модели
minalt = models[modelI].altmin
maxalt = models[modelI].altmax

# Если необходимо, измените диапазоны в соответствии с координатами
if igdgc == 2:
    minalt += 6371.2
    maxalt += 6371.2

# Вызов функции для расчета магнитных компонентов
calculate(sdate, modelI, igdgc, latitude, longitude, alt, x, y, z, d, f, h, i, ddot, fdot, hdot, idot, xdot, ydot, zdot)

# Вывод результатов
zdot = 1


import math

def calculate(sdate, modelI, igdgc, latitude, longitude, alt, x, y, z, d, f, h, i, ddot, fdot, hdot, idot, xdot, ydot, zdot):
    print(f"{sdate:.5f} {latitude:.5f} {longitude:.5f} {alt:.5f}")

    dtemp, ftemp, htemp, itemp = 0, 0, 0, 0
    xtemp, ytemp, ztemp = 0, 0, 0

    if max2[modelI] == 0:
        getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1)
        getshc(mdfile, 1, irec_pos[modelI+1], max1[modelI+1], 2)
        nmax = interpsh(sdate, yrmin[modelI], max1[modelI],
                        yrmin[modelI+1], max1[modelI+1], 3)
        nmax = interpsh(sdate+1, yrmin[modelI], max1[modelI],
                        yrmin[modelI+1], max1[modelI+1], 4)
    else:
        getshc(mdfile, 1, irec_pos[modelI], max1[modelI], 1)
        getshc(mdfile, 0, irec_pos[modelI], max2[modelI], 2)
        nmax = extrapsh(sdate, epoch[modelI], max1[modelI], max2[modelI], 3)
        nmax = extrapsh(sdate+1, epoch[modelI], max1[modelI], max2[modelI], 4)

    shval3(igdgc, latitude, longitude, alt, nmax, 3,
           IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3,
           x, y, z, xtemp, ytemp, ztemp)

    dihf(3, x, y, z, xtemp, ytemp, ztemp, d, f, h, i, dtemp, ftemp, htemp, itemp)

    shval3(igdgc, latitude, longitude, alt, nmax, 4,
           IEXT, EXT_COEFF1, EXT_COEFF2, EXT_COEFF3,
           x, y, z, xtemp, ytemp, ztemp)

    dihf(4, x, y, z, xtemp, ytemp, ztemp, d, f, h, i, dtemp, ftemp, htemp, itemp)

    ddot = (dtemp - d) * math.degrees(1)
    if ddot > 180.0:
        ddot -= 360.0
    if ddot <= -180.0:
        ddot += 360.0
    ddot *= 60.0

    idot = (itemp - i) * math.degrees(1) * 60
    d = d * math.degrees(1)
    i = i * math.degrees(1)
    hdot = htemp - h
    xdot = xtemp - x
    ydot = ytemp - y
    zdot = ztemp - z
    fdot = ftemp - f

    if h < 100.0:
        d = float("nan")
        ddot = float("nan")

    if h < 1000.0:
        warn_H = 0
        warn_H_strong = 1
        if h < warn_H_strong_val:
            warn_H_strong_val = h
    elif h < 5000.0 and not warn_H_strong:
        warn_H = 1
        if h < warn_H_val:
            warn_H_val = h

    if abs(90.0 - abs(latitude)) <= 0.001:
        x = float("nan")
        y = float("nan")
        d = float("nan")
        xdot = float("nan")
        ydot = float("nan")
        ddot = float("nan")
        warn_P = 1
        warn_H = 0
        warn_H_strong = 0
