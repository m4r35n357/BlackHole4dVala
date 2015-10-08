#!/usr/bin/python

from __future__ import absolute_import, division, print_function, unicode_literals
from numpy import sqrt, sin, cos, radians, pi, add, subtract, array
import demo
import pi3d
from sys import argv, stdout, stderr
from json import loads
from pprint import pprint

class Planet(pi3d.Sphere):
  def __init__(self, shader, colour, radius, pos=[0.0, 0.0, 0.0], track_shader=None):
    super(Planet, self).__init__(radius=radius, x=pos[0], y=pos[1], z=pos[2])
    super(Planet, self).set_draw_details(shader, [])
    self.pos = array(pos)
    self.track_shader = track_shader
    self.t_v = []
    self.t_len = 0
    self.trace_shape = None
    self.set_material(colour)
    
  def position_and_draw(self):
    self.position(self.pos[0], self.pos[1], self.pos[2])
    self.draw()
    if self.track_shader != None:
      self.t_len += 1
      self.t_v.append(tuple(self.pos))
      if (self.t_len % 10) == 0:
        self.trace_shape = pi3d.Points(vertices=self.t_v, material=(1.0,0.0,0.0), point_size=1)
        self.trace_shape.set_shader(self.track_shader)
        if (self.t_len > 2400):
          self.t_v = self.t_v[-2400:]
          self.t_len = 2400
      if self.trace_shape:
        self.trace_shape.draw()

if len(argv) < 3:
    raise Exception('>>> ERROR! Please supply a data file name and a parameter file name <<<')
dataFile = open(argv[1], 'r')
line = dataFile.readline()
parameters = loads(open(argv[2], 'r').read())
a = parameters['a']
a2 = a**2
m = parameters['M']
horizon = m * (1.0 + sqrt(1.0 - a2))
# Setup display and initialise pi3d ------
DISPLAY = pi3d.Display.create(x=0, y=0, frames_per_second=0)
DISPLAY.set_background(0,0,0,1)    	# r,g,b,alpha
# Camera ---------------------------------
CAMERA = pi3d.Camera()
rot = 0
tilt = 0
rottilt = True
camRad = 30.0
# Shaders --------------------------------
shader = pi3d.Shader("mat_light")
tracksh = pi3d.Shader("mat_flat")
# Planets --------------------------------
sun = Planet(shader, (0.0,0.0,1.0), horizon, pos=[0.0, 0.0, 0.0])
earth = Planet(shader, (0.0,1.0,0.0), 0.125, track_shader=tracksh)
# Fetch key presses ----------------------
mykeys = pi3d.Keyboard()
# Display scene
while DISPLAY.loop_running():
  params = loads(line)
  if rottilt:
    CAMERA.reset()
    CAMERA.rotate(-tilt, rot, 0)
    CAMERA.position((camRad * sin(radians(rot)) * cos(radians(tilt)), camRad * sin(radians(tilt)), -camRad * cos(radians(rot)) * cos(radians(tilt))))
    rottilt = False
  sun.position_and_draw()
  r = float(params['r'])
  th = float(params['th'])
  ph = float(params['ph'])
  rasth = sqrt(r**2 + a2) * sin(th)
  earth.pos = [rasth * cos(ph), rasth * sin(ph), r * cos(th)]
  earth.position_and_draw()
  k = mykeys.read()
  if k >-1:
    rottilt = True
    if k==112:
      pi3d.screenshot("orbit.jpg")
    elif k==119:  #key W rotate camera up
      tilt += 2.0
    elif k==115:  #kry S down
      tilt -= 2.0
    elif k==97:   #key A left
      rot -= 2
    elif k==100:  #key D right
      rot += 2
    elif k==61:   #key += in
      camRad -= 0.5
    elif k==45:   #key _- out
      camRad += 0.5
    elif k==27:
      mykeys.close()
      DISPLAY.destroy()
      break
  line = dataFile.readline()
  if not line:
    DISPLAY.stop()

