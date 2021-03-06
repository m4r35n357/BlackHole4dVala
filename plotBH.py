#!/usr/bin/env python
"""
Copyright (c) 2014, 2015, 2016, 2017 Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from json import loads
from math import sqrt, sin, cos, fabs, pi, atan2, acos, log10
from sys import argv, stderr, stdin
from visual import display, sphere, curve, rate, ellipsoid, ring, color, label, points


def isco (a, l):
    z1 = 1.0 + pow(1.0 - a**2, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0))
    z2 = sqrt(3.0 * a**2 + z1 * z1)
    if a * l >= 0.0:
        return 3.0 + z2 - sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))
    else:
        return 3.0 + z2 + sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))

def p_sphere (a, l):
    if a * l >= 0.0:
        return 2.0 * (1.0 + cos(2.0 / 3.0 * acos(-a if a >= 0.0 else a)))
    else:
        return 2.0 * (1.0 + cos(2.0 / 3.0 * acos(a if a >= 0.0 else -a)))

def speed (gamma):
    if gamma > 1.0 or -gamma > 1.0:
        return sqrt(1.0 - 1.0 / gamma**2)
    else:
        return 0.0
    #return sqrt(gamma * (1.0 - 2.0 * r / (r**2 + a**2 * cos(theta)**2)) - 1.0 / gamma**2)

def to_rectangular (polar, a):
    ra = sqrt(polar[0]**2 + a**2)
    sth = sin(polar[1])
    return ra * sth * cos(polar[2]), ra * sth * sin(polar[2]), polar[0] * cos(polar[1])

def error_colour (error):
    if error < 1.0e-18:
        return color.blue
    elif error < 1.0e-12:
        return color.green
    elif error < 1.0e-9:
        return color.cyan
    elif error < 1.0e-6:
        return color.yellow
    elif error < 1.0e-3:
        return color.orange
    else:
        return color.red

def log_error(e):
    """
    convert error to a pseudo db Value
    :param e: the numerical value of he error
    :return: the dB value, clamped to a minimum
    """
    error = e if e >= 0.0 else -e
    return 10.0 * log10(error) if error > 1.0e-36 else -360.0


def main():
    print("Geodesic Plotter: {}".format(argv))
    if len(argv) < 2:
        raise Exception('>>> ERROR! Please supply a parameter file name <<<')
    parameters = loads(open(argv[1]).read())['IC']
    m = parameters['M'] if 'M' in parameters else 1.0
    a = parameters['a'] if 'a' in parameters else  1.0
    l = parameters['L'] if 'L' in parameters else 0.0
    mu = parameters['mu'] if 'mu' in parameters else 0.0
    horizon = m * (1.0 + sqrt(1.0 - a**2))
    cauchy = m * (1.0 - sqrt(1.0 - a**2))
    #  set up the scene
    window_name = 'orbit'
    my_scene = display(title=window_name)
    my_scene.center = centre = (0.0, 0.0, 0.0)
    my_scene.width = my_scene.height = 1024
    my_scene.range = (20.0, 20.0, 20.0)
    # Inner Horizon
    inner = 2.0 * to_rectangular((cauchy, 0.5 * pi, 0.0), a)[0]
    ellipsoid(pos=centre, length=inner, height=inner, width=2.0*to_rectangular((cauchy, 0.0, 0.0), a)[2],
              color=color.blue, opacity=0.6)
    # Outer Horizon
    outer = 2.0 * to_rectangular((horizon, 0.5 * pi, 0.0), a)[0]
    ellipsoid(pos=centre, length=outer, height=outer, width=2.0*to_rectangular((horizon, 0.0, 0.0), a)[2],
              color=color.blue, opacity=0.4)
    # Ergosphere
    ergo = 2.0 * to_rectangular((2.0 * m, 0.5 * pi, 0.0), a)[0]
    ellipsoid(pos=centre, length=ergo, height=ergo, width=2.0*to_rectangular((horizon, 0.0, 0.0), a)[2],
              color=color.gray(0.7), opacity=0.2)
    # Singularity
    if fabs(a) > 0.0:
        ring(pos=centre, axis=(0, 0, 1), radius=a, color=color.white, thickness=0.01)
    else:
        sphere(pos=centre, radius=0.05, color=color.white)
    # ISCO
    ring(pos=centre, axis=(0, 0, 1), radius=to_rectangular((isco(a, l), 0.5 * pi, 0.0), a)[0],
         color=color.magenta, thickness=0.01)
    # Photon Sphere(s)
    ring(pos=centre, axis=(0, 0, 1), radius=to_rectangular((p_sphere(a, l), 0.5 * pi, 0.0), a)[0],
         color=color.orange, thickness=0.01)
    ring(pos=centre, axis=(0, 0, 1), radius=to_rectangular((p_sphere(-a, l), 0.5 * pi, 0.0), a)[0],
         color=color.orange, thickness=0.01)
    #photons = 2.0 * to_rectangular((pSphere(a, l), 0.5 * pi, 0.0), a)[0]
    #ellipsoid(pos=centre, length=photons, height=photons, width=2.0*to_rectangular((pSphere(a, l), 0.0, 0.0), a)[2],
    #     color=color.white, opacity=0.1)
    # z axis
    curve(pos=[(0.0, 0.0, -15.0), (0.0, 0.0, 15.0)], color = color.gray(0.7))
    # radial line
    radial = curve(pos=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], color=color.gray(0.2))
    # Data display
    hud = label(pos=(0.0, 0.0, 0.0), xoffset=340, yoffset=330, line=False, border=10, font='Monospace', height=16,
                    color=(0.5, 0.5, 0.0), linecolor=color.gray(0.1))
    #cone(pos=(0,0,12), axis=(0,0,-12), radius=12.0 * tan(0.15 * pi), opacity=0.2)
    #cone(pos=(0,0,-12), axis=(0,0,12), radius=12.0 * tan(0.15 * pi), opacity=0.2)
    #sphere(pos=(0,0,0), radius=3.0, opacity=0.2)
    #sphere(pos=(0,0,0), radius=12.0, opacity=0.1)
    # animate!
    ball = sphere()  # Particle
    count = 0
    data_line = stdin.readline()
    t_old = 0.0
    e_cum = e_pk = 0.0
    while data_line:  # build raw data arrays
        rate(60)
        if count % 1000 == 0:
            ball.visible = False
            ball = sphere(radius=0.2)  # Particle
            ball.trail = curve(size=1)  #  trail
        data = loads(data_line)
        error = data['v4e']
        e = error if error >= 0.0 else -error
        count += 1
        ball.color = error_colour(e_cum / count)
        r = data['r']
        th = data['th']
        ph = data['ph']
        ball.pos = to_rectangular((r, th, ph), a)
        radial.pos = ((0.0, 0.0, 0.0), ball.pos)
        if data['tP'] * t_old < 0 or data['tP'] - t_old > 100.0:
            ball.color = color.white
        ball.trail.append(pos=ball.pos, color=error_colour(e))
        # Data display
        if fabs(mu) > 0.0:
            hud.text = u"v  %.6f\n\u03c4  %.1f\nt  %.1f\nr  %.3f\n\u03b8  %.0f\n\u03d5  %.0f" % (speed(data['tP']), data['tau'], data['t'],
                                                                                                 r, atan2(ball.pos.z, sqrt(ball.pos.x**2 + ball.pos.y**2)) * 180.0 / pi, ph * 180.0 / pi % 360.0)
        else:
            hud.text = u"\u03bb  %.1f\nt  %.1f\nr  %.3f\n\u03b8  %.0f\n\u03d5  %.0f" % (data['tau'], data['t'],
                                                                                        r, atan2(ball.pos.z, sqrt(ball.pos.x**2 + ball.pos.y**2)) * 180.0 / pi, ph * 180.0 / pi % 360.0)
        #popen('import -window ' + windowName + ' -compress None VPythonOutput/' + str(counter).zfill(4) + '.png')
        t_old = data['tP']
        e_cum += e
        e_pk = e_pk if e_pk > e else e
        data_line = stdin.readline()
    print(" {}: Errors - Peak: {:.1f}, Average: {:.1f}".format(argv[0], log_error(e_pk), log_error(e_cum / count)))

if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
