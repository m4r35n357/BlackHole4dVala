#!/usr/bin/python3

from __future__ import absolute_import, division, print_function, unicode_literals

from json import loads
from numpy import sqrt, sin, cos, radians, array
from sys import argv, stdin, stderr

from pi3d import Sphere, Points, Display, Camera, Shader, Keyboard, screenshot, Lines


def error_colour (error):
    e = abs(error)
    if e < 1.0e-12:
        return 0.0, 1.0, 0.0
    elif e < 1.0e-9:
        return 0.0, 0.5, 0.5
    elif e < 1.0e-6:
        return 1.0, 1.0, 0.0
    elif e < 1.0e-3:
        return 0.5, 0.5, 0.0
    else:
        return 1.0, 0.0, 0.0

class Body(Sphere):
    def __init__(self, shader, colour, radius, position=(0.0, 0.0, 0.0), track_shader=None):
        super(Body, self).__init__(radius=radius, x=position[0], y=position[1], z=position[2])
        super(Body, self).set_draw_details(shader, [])
        self.pos = array(position)
        self.track_shader = track_shader
        self.vertices = []
        self.track_length = 0
        self.track_max = 2400
        self.trace_shape = None
        self.set_material(colour)

    def position_and_draw(self, trace_material=(0.5, 0.5, 0.5)):
        # body
        self.position(self.pos[0], self.pos[1], self.pos[2])
        self.draw()
        # track
        if self.track_shader:
            self.track_length += 1
            self.vertices.append(tuple(self.pos))
            if (self.track_length % 1) == 0:
                self.trace_shape = Lines(vertices=self.vertices, material=trace_material)
                self.trace_shape.set_shader(self.track_shader)
                if self.track_length > self.track_max:
                    self.vertices = self.vertices[-self.track_max:]
                    self.track_length = self.track_max
            if self.trace_shape:
                self.trace_shape.draw()

def main():
    print("pi3d Geodesic Plotter: {}".format(argv))
    if len(argv) < 2:
        raise Exception('>>> ERROR! Please supply a parameter file name <<<')
    parameters = loads(open(argv[1]).read())['IC']
    interval = parameters['plotratio']
    m = parameters['M']
    a = parameters['a']
    a2 = a**2
    # Setup display and initialise pi3d
    display = Display.create(x=0, y=0, frames_per_second=0)
    display.set_background(0, 0, 0, 1)  # r,g,b,alpha
    # Camera
    camera = Camera()
    rot = tilt = 0
    rot_tilt = True
    cam_rad = 30.0
    # Bodies
    body_shader = Shader("mat_light")
    track_shader = Shader("mat_flat")
    black_hole = Body(body_shader, (0.0, 0.0, 1.0), m * (1.0 + sqrt(1.0 - a2)), position=[0.0, 0.0, 0.0])
    particle = Body(body_shader, (0.0, 1.0, 0.0), 0.125, track_shader=track_shader)
    # Enable key presses
    axis = Lines(vertices=[(0, 0, 5,), (0, 0, -5,)], line_width=10)
    axis.set_draw_details(track_shader, [])
    keys = Keyboard()
    # Display scene
    counter = 1
    cumulative_error = 0.0
    line = stdin.readline()
    while display.loop_running():
        data = loads(line)
        # monitor errors
        current_error = data['v4e']
        cumulative_error += current_error if current_error >= 0.0 else -current_error
        # camera control
        if rot_tilt:
            camera.reset()
            camera.rotate(-tilt, rot, 0)
            camera.position((cam_rad * sin(radians(rot)) * cos(radians(tilt)), cam_rad * sin(radians(tilt)),
                             -cam_rad * cos(radians(rot)) * cos(radians(tilt))))
            rot_tilt = False
        # plot the black hole
        black_hole.position_and_draw()
        axis.draw()
        # plot the orbiter
        if counter % interval == 0:
            r = float(data['r'])
            th = float(data['th'])
            ph = float(data['ph'])
            ra_sth = sqrt(r**2 + a2) * sin(th)
            particle.pos = [ra_sth * cos(ph), ra_sth * sin(ph), r * cos(th)]
            particle.set_material(error_colour(data['v4e']))
            particle.position_and_draw(trace_material=error_colour(cumulative_error / counter))
        # process keyboard input
        key = keys.read()
        if key > -1:
            rot_tilt = True
            if key == 112:
                screenshot("orbit.jpg")
            elif key == 119:  # key W rotate camera up
                tilt += 2.0
            elif key == 115:  # kry S down
                tilt -= 2.0
            elif key == 97:  # key A left
                rot -= 2
            elif key == 100:  # key D right
                rot += 2
            elif key == 61:  # key += in
                cam_rad -= 0.5
            elif key == 45:  # key _- out
                cam_rad += 0.5
            elif key == 27:
                keys.close()
                display.destroy()
                break
        # prepare for next iteration
        line = stdin.readline()
        counter += 1
        if not line:
            display.stop()

if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
