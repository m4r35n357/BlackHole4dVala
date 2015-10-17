#!/usr/bin/python

from __future__ import absolute_import, division, print_function, unicode_literals
from numpy import sqrt, sin, cos, radians, pi, add, subtract, array
import pi3d
from sys import argv, stdout, stderr
from json import loads
from pprint import pprint

def isco (a):
	z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0))
	z2 = sqrt(3.0 * a * a + z1 * z1)
	if a >= 0.0:
		return 3.0 + z2 - sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))
	else:
		return 3.0 + z2 + sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))

class Planet(pi3d.Sphere):
  def __init__(self, shader, colour, radius, pos=[0.0, 0.0, 0.0], track_shader=None, scale=1.0, grid=12, opacity=1.0):
    super(Planet, self).__init__(radius=radius, slices=grid, sides=grid, sx=scale, sy=scale, sz=1.0, x=pos[0], y=pos[1], z=pos[2])
    self.set_draw_details(shader, [])
    self.pos = array(pos)
    self.track_shader = track_shader
    self.t_v = []
    self.t_len = 0
    self.trace_shape = None
    self.set_material(colour)
    self.set_alpha(opacity)
    
  def position_and_draw(self, trace_material=(0.5,0.5,0.5)):
    self.position(self.pos[0], self.pos[1], self.pos[2])
    self.draw()
    if self.track_shader != None:
      self.t_len += 1
      self.t_v.append(tuple(self.pos))
      if (self.t_len % 10) == 0:
        self.trace_shape = pi3d.Points(vertices=self.t_v, material=trace_material, point_size=1)
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
cauchy = m * (1.0 - sqrt(1.0 - a2))
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
bulge = sqrt(1.0 + a2)
ergo = Planet(shader, (1.0,1.0,0.0), horizon, pos=[0.0, 0.0, 0.0], scale=sqrt(2.0 + a2), grid=24, opacity=0.1)
outer = Planet(shader, (0.0,1.0,1.0), horizon, pos=[0.0, 0.0, 0.0], scale=bulge, grid=24, opacity=0.3)
inner = Planet(shader, (0.0,0.0,1.0), cauchy, pos=[0.0, 0.0, 0.0], scale=bulge, grid=24)
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
  inner.position_and_draw()
  outer.position_and_draw()
  ergo.position_and_draw()
  r = float(params['r'])
  th = float(params['th'])
  ph = float(params['ph'])
  rasth = sqrt(r**2 + a2) * sin(th)
  earth.pos = [rasth * cos(ph), rasth * sin(ph), r * cos(th)]
  error = params['v4e']
  if error < -120.0:
    colour = (0.0,1.0,0.0)
  elif error < -90.0:
    colour = (0.0,0.5,0.5)
  elif error < -60.0:
    colour = (1.0,1.0,0.0)
  elif error < -30.0:
    colour = (0.5,0.5,0.0)
  else:
    colour = (1.0,0.0,0.0)
  earth.set_material(colour)
  earth.position_and_draw(trace_material=colour)
#  Cylinder(radius=isco(a), height=0.1).draw()
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

