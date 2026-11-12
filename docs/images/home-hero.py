"""Draws the figure of the home page, docs/images/home-hero.webp, from a solution of the HyShot II case.

    python docs/images/home-hero.py <field.plt> [out.webp]

<field.plt> is the field file that MOSE writes for test/3D/hyshotII (Tecplot binary, cell-centred). The computed
half pitch is mirrored into four injectors; the figure shows the hydrogen of the jets (Y_H2 = 0.08, cyan), the flame
around them (Y_OH = 0.005, coloured by temperature, 1300 to 2500 K) and the bottom wall (coloured by pressure), fading
out upstream of the injectors and towards x = 660 mm. Needs numpy, matplotlib, Pillow and PyVista, and an OpenGL
context: without DISPLAY the script starts its own Xvfb, since a VTK with neither EGL nor OSMesa crashes without one.
"""
import os
import struct
import subprocess
import sys

import matplotlib
import matplotlib.colors as mc
import numpy as np
import pyvista as pv
from PIL import Image

W = 0.009375                 # half pitch of the injectors [m]
NINJ = 4                     # injectors drawn
XMIN, X0, X1 = 0.396, 0.61, 0.71
UFADE = 0.012                # length of the fade upstream [m]
WALL_MAP = mc.LinearSegmentedColormap.from_list('wall', ['#06111f', '#0e2a4a', '#1f5f9e', '#5fb4f0', '#d8efff'])


def read_plt(path, want):
    """The ordered zones of a Tecplot binary file (#!TDV112): per zone, the variables of `want` as (K, J, I) arrays,
    cell-centred ones as (K-1, J-1, I-1)."""
    buf = memoryview(np.memmap(path, dtype=np.uint8, mode='r'))
    p = 12

    def i32():
        nonlocal p
        p += 4
        return struct.unpack_from('<i', buf, p - 4)[0]

    def f32():
        nonlocal p
        p += 4
        return struct.unpack_from('<f', buf, p - 4)[0]

    def text():
        s = []
        while (c := i32()) != 0:
            s.append(chr(c))
        return ''.join(s)

    i32(); text()
    nv = i32()
    names = [text() for _ in range(nv)]
    zones = []
    while abs((m := f32()) - 357.0) > 1e-3:
        assert abs(m - 299.0) < 1e-3, 'unexpected record in the header'
        text(); i32(); i32(); p += 12
        assert i32() == 0, 'not an ordered zone'
        loc = [i32() for _ in range(nv)] if i32() == 1 else [0] * nv
        i32(); assert i32() == 0
        dims = (i32(), i32(), i32())
        assert i32() == 0
        zones.append((dims, loc))
    out = []
    for (I, J, K), loc in zones:
        assert abs(f32() - 299.0) < 1e-3
        fmt = [i32() for _ in range(nv)]
        assert i32() == 0 and i32() == 0, 'passive or shared variables'
        i32()
        p += 16 * nv
        data = {}
        for v in range(nv):
            n = I * J * K if loc[v] == 0 else I * J * (K - 1)     # cell-centred values are padded in I and J
            dt = np.dtype({1: '<f4', 2: '<f8'}[fmt[v]])
            if names[v] in want:
                a = np.frombuffer(buf, dtype=dt, count=n, offset=p)
                data[names[v]] = a.reshape(K, J, I) if loc[v] == 0 else a.reshape(K - 1, J, I)[:, :J - 1, :I - 1]
            p += n * dt.itemsize
        out.append(data)
    return out


def surfaces(path):
    """Iso-surfaces of the jets and of the flame, and the bottom wall, of the computed half pitch."""
    want = {'x', 'y', 'z', 'T', 'p'} | {f'rho({i})' for i in range(1, 10)}
    h2, oh, wall = [], [], []
    for d in read_plt(path, want):
        g = pv.StructuredGrid(*(d[c].transpose(2, 1, 0) for c in 'xyz'))
        rho = sum(d[f'rho({i})'] for i in range(1, 10))
        g.cell_data['T'] = d['T'].ravel()
        g.cell_data['p'] = d['p'].ravel()
        g.cell_data['YH2'] = (d['rho(3)'] / rho).ravel()   # species: N2 O2 H2 H2O OH O H HO2 H2O2
        g.cell_data['YOH'] = (d['rho(5)'] / rho).ravel()
        gp = g.cell_data_to_point_data()
        for f, level, lst in (('YH2', 0.45, h2), ('YOH', 0.01, oh)):
            c = gp.contour([level], scalars=f)
            if c.n_points:
                lst.append(c)
        s = g.extract_surface(algorithm='dataset_surface')
        normal = s.compute_normals(cell_normals=True, point_normals=False).cell_data['Normals'][:, 2]
        bottom = np.where((s.cell_centers().points[:, 2] < 1e-6) & (np.abs(normal) > 0.99))[0]
        if bottom.size:
            wall.append(s.extract_cells(bottom).extract_surface(algorithm='dataset_surface').cell_data_to_point_data())
    return [lst[0].merge(lst[1:]).extract_surface(algorithm='dataset_surface') for lst in (h2, oh, wall)]


def replicate(m):
    """The half pitch mirrored into NINJ injectors, clipped and given a fade (0 to 1) along x."""
    r = m.reflect((0, 1, 0), point=(0, 0, 0))
    parts = []
    for k in range(NINJ):
        dy = (k - (NINJ - 1) / 2) * 2 * W
        parts += [m.translate((0, dy, 0), inplace=False), r.translate((0, dy, 0), inplace=False)]
    out = parts[0].merge(parts[1:]).clip_box((XMIN, X1, -1, 1, -1, 1), invert=False).extract_surface(algorithm='dataset_surface')
    x = out.points[:, 0]
    s = np.clip((X1 - x) / (X1 - X0), 0, 1) * np.clip((x - XMIN) / UFADE, 0, 1)
    out.point_data['fade'] = s * s * (3 - 2 * s)
    return out


def paint(m, scalars, cmap, clim, opacity):
    """Colours and opacity per point, so that the fade applies to the coloured surfaces too."""
    if scalars is None:
        rgb = np.tile(mc.to_rgb(cmap), (m.n_points, 1))
    else:
        cmap = matplotlib.colormaps[cmap] if isinstance(cmap, str) else cmap
        rgb = cmap(np.clip((m[scalars] - clim[0]) / (clim[1] - clim[0]), 0, 1))[:, :3]
    m.point_data['rgba'] = (np.column_stack([rgb, m['fade'] * opacity]) * 255).astype(np.uint8)
    return m


def start_xvfb():
    """A private Xvfb, on the display it picks itself, when DISPLAY is not set; None otherwise."""
    if os.environ.get('DISPLAY'):
        return None
    r, w = os.pipe()
    try:
        xvfb = subprocess.Popen(['Xvfb', '-displayfd', str(w), '-screen', '0', '1920x1200x24', '+extension', 'GLX',
                                 '-nolisten', 'tcp'], pass_fds=(w,), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except FileNotFoundError:
        sys.exit('home-hero.py: DISPLAY is not set and Xvfb is not installed: VTK needs a display to render')
    os.close(w)
    with os.fdopen(r) as f:
        number = f.readline().strip()        # written by Xvfb once it accepts connections
    if not number:
        sys.exit('home-hero.py: Xvfb did not start')
    os.environ['DISPLAY'] = ':' + number
    return xvfb


def main(path, out='docs/images/home-hero.webp'):
    xvfb = start_xvfb()
    try:
        render(path, out)
    finally:
        if xvfb:
            xvfb.terminate()


def render(path, out):
    h2, oh, wall = (replicate(m) for m in surfaces(path))
    p = pv.Plotter(off_screen=True, window_size=(2800, 1800))
    p.set_background('black')        # VTK blends with the background: on black the pixels come premultiplied (RGBa)
    p.add_mesh(paint(wall, 'p', WALL_MAP, (4e4, 3e5), 1.0), scalars='rgba', rgba=True,
               ambient=0.3, diffuse=0.7, specular=0.1)
    p.add_mesh(paint(h2, None, '#9fe3ff', None, 1.0), scalars='rgba', rgba=True,
               smooth_shading=True, specular=0.6, specular_power=30)
    p.add_mesh(paint(oh, 'T', 'turbo', (1500, 2200), 0.55), scalars='rgba', rgba=True,
               smooth_shading=True, ambient=0.25, specular=0.3)
    p.enable_anti_aliasing('ssaa')
    cam = p.camera
    cam.focal_point, cam.position, cam.up, cam.view_angle = (0.49, 0, 0), (0.298, -0.142, 0.132), (0, 0, 1), 30
    shot = np.ascontiguousarray(p.screenshot(None, transparent_background=True, return_img=True))
    p.close()                        # releases the X connection before a private Xvfb stops
    img = Image.frombytes('RGBa', (shot.shape[1], shot.shape[0]), shot.tobytes()).convert('RGBA')
    alpha = np.asarray(img)[:, :, 3]
    ys, xs = np.nonzero(alpha > 4)
    img = img.crop((max(xs.min() - 24, 0), max(ys.min() - 24, 0), min(xs.max() + 25, img.width), min(ys.max() + 25, img.height)))
    img.thumbnail((1600, 1600), Image.LANCZOS)
    img.save(out, quality=90, method=6)
    print(out, img.size)


if __name__ == '__main__':
    main(*sys.argv[1:])
