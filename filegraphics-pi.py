#!/usr/bin/python3

from __future__ import absolute_import, division, print_function, unicode_literals

from json import loads
from numpy import sqrt, sin, cos, radians, array
from sys import argv, stdin

from pi3d import Sphere, Points, Display, Camera, Shader, Keyboard, screenshot


def error_colour (e):
    error = abs(e)
    if error < 1.0e-12:
        return 0.0, 1.0, 0.0
    elif error < 1.0e-9:
        return 0.0, 0.5, 0.5
    elif error < 1.0e-6:
        return 1.0, 1.0, 0.0
    elif error < 1.0e-3:
        return 0.5, 0.5, 0.0
    else:
        return 1.0, 0.0, 0.0


class Planet(Sphere):
    def __init__(self, shdr, clr, radius, pos=(0.0, 0.0, 0.0), track_shader=None):
        super(Planet, self).__init__(radius=radius, x=pos[0], y=pos[1], z=pos[2])
        super(Planet, self).set_draw_details(shdr, [])
        self.pos = array(pos)
        self.track_shader = track_shader
        self.t_v = []
        self.t_len = 0
        self.trace_shape = None
        self.set_material(clr)

    def position_and_draw(self, trace_material=(0.5, 0.5, 0.5)):
        self.position(self.pos[0], self.pos[1], self.pos[2])
        self.draw()
        if self.track_shader is not None:
            self.t_len += 1
            self.t_v.append(tuple(self.pos))
            if (self.t_len % 10) == 0:
                self.trace_shape = Points(vertices=self.t_v, material=trace_material)
                self.trace_shape.set_shader(self.track_shader)
                if self.t_len > 2400:
                    self.t_v = self.t_v[-2400:]
                    self.t_len = 2400
            if self.trace_shape:
                self.trace_shape.draw()


if len(argv) < 2:
    raise Exception('>>> ERROR! Please supply a parameter file name <<<')
parameters = loads(open(argv[1]).read())
interval = parameters['IC']['plotratio']
m = parameters['IC']['M']
a = parameters['IC']['a']
a2 = a**2
horizon = m * (1.0 + sqrt(1.0 - a2))
# Setup display and initialise pi3d ------
DISPLAY = Display.create(x=0, y=0, frames_per_second=0)
DISPLAY.set_background(0, 0, 0, 1)  # r,g,b,alpha
# Camera ---------------------------------
CAMERA = Camera()
rot = tilt = 0
rot_tilt = True
camRad = 30.0
# Shaders --------------------------------
shader = Shader("mat_light")
tracksh = Shader("mat_flat")
# Planets --------------------------------
sun = Planet(shader, (0.0, 0.0, 1.0), horizon, pos=[0.0, 0.0, 0.0])
earth = Planet(shader, (0.0, 1.0, 0.0), 0.125, track_shader=tracksh)
# Fetch key presses ----------------------
mykeys = Keyboard()
# Display scene
colour = None
counter = 0
line = stdin.readline()
while DISPLAY.loop_running():
    params = loads(line)
    if rot_tilt:
        CAMERA.reset()
        CAMERA.rotate(-tilt, rot, 0)
        CAMERA.position((camRad * sin(radians(rot)) * cos(radians(tilt)), camRad * sin(radians(tilt)),
                         -camRad * cos(radians(rot)) * cos(radians(tilt))))
        rot_tilt = False
    sun.position_and_draw()
    if counter % interval == 0:
        r = float(params['r'])
        th = float(params['th'])
        ph = float(params['ph'])
        ra_sth = sqrt(r**2 + a2) * sin(th)
        earth.pos = [ra_sth * cos(ph), ra_sth * sin(ph), r * cos(th)]
        colour = error_colour( params['v4e'])
        earth.set_material(colour)
    earth.position_and_draw(trace_material=colour)
    k = mykeys.read()
    if k > -1:
        rot_tilt = True
        if k == 112:
            screenshot("orbit.jpg")
        elif k == 119:  # key W rotate camera up
            tilt += 2.0
        elif k == 115:  # kry S down
            tilt -= 2.0
        elif k == 97:  # key A left
            rot -= 2
        elif k == 100:  # key D right
            rot += 2
        elif k == 61:  # key += in
            camRad -= 0.5
        elif k == 45:  # key _- out
            camRad += 0.5
        elif k == 27:
            mykeys.close()
            DISPLAY.destroy()
            break
    line = stdin.readline()
    counter += 1
    if not line:
        DISPLAY.stop()
