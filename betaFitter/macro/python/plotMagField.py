#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import uproot
from mayavi import mlab


def read_magnetic_field(root_file):
    with uproot.open(root_file) as file:
        data = file["magfield"].arrays(["x", "y", "z", "bx", "by", "bz", "b_magnitude"], library="np")
        return (data[f].reshape(150,150,150) for f in ["x", "y", "z", "bx", "by", "bz", "b_magnitude"])

x, y, z, u, v, w, mag = read_magnetic_field("mag_field.root")

src = mlab.pipeline.vector_field(x, y, z, u, v, w, scalars=mag)

streamlines = mlab.pipeline.streamline(
    src, 
    seedtype='plane',
    seed_visible=True,
    seed_resolution=20,
    integration_direction='both'
)

streamlines.seed.widget.origin = [-150, -150, 0]
streamlines.seed.widget.point1 = [150, -150, 0]
streamlines.seed.widget.point2 = [-150, 150, 0]

streamlines.stream_tracer.maximum_propagation = 500
streamlines.tube_filter.radius = 0.4
streamlines.streamline_type = 'tube'

mlab.colorbar(streamlines, title='Magnetic Field', orientation='vertical')
axes = mlab.axes(xlabel='X', ylabel='Y', zlabel='Z', nb_labels=5)
mlab.title('Uniform Magnetic Field Streamlines')

mlab.show()