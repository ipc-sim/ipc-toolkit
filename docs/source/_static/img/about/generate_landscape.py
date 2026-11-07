"""Generate the "landscape" figure shown on about/why.rst.

Run this script to rewrite ``_static/landscape.html`` after editing POINTS:

    python docs/source/_static/img/about/generate_landscape.py

We emit inline SVG (embedded with ``.. raw:: html``) instead of a standalone
``.svg`` file. A standalone file loaded through ``<img>`` cannot see the theme's
CSS variables, so it would not follow the light/dark toggle. Inline SVG picks
up the colors from ``custom.css`` and switches with the rest of the page.

The positions are illustrative. They come from each project's stated scope, not
from a benchmark, so nudge a point if you think it belongs somewhere else.
"""

import pathlib
from xml.sax.saxutils import escape

# Both axes run from 0 to 1.
#   x: 0 = modular component you embed, 1 = integrated framework you build on
#   y: 0 = generic contact handling,    1 = IPC-grade (barrier + CCD)
#
# kind is one of:
#   "toolkit" -- the IPC Toolkit itself
#   "uses"    -- projects built on top of the IPC Toolkit
#   "other"   -- everything else
#
# anchor says which side of the dot the label goes on: "start" puts the label
# to the right of the dot, "end" puts it to the left. Labels that would collide
# with a neighbor get a dy offset (in pixels) to slide them up or down.
POINTS = [
    # Modular components
    dict(name="IPC Toolkit", x=0.13, y=0.86, kind="toolkit", anchor="start"),
    dict(name="Tight-Inclusion / Scalable-CCD", x=0.07, y=0.60, kind="other", anchor="start"),
    dict(name="FCL / Coal", x=0.17, y=0.42, kind="other", anchor="start"),
    dict(name="Ipopt", x=0.30, y=0.27, kind="other", anchor="start"),
    dict(name="TinyAD / Bartels", x=0.16, y=0.14, kind="other", anchor="start"),
    dict(name="SimKit", x=0.42, y=0.11, kind="other", anchor="start"),
    # Integrated frameworks with IPC-grade contact
    dict(name="libuipc", x=0.79, y=0.93, kind="other", anchor="end"),
    dict(name="GIPC", x=0.92, y=0.81, kind="other", anchor="end"),
    dict(name="IPC / C-IPC", x=0.95, y=0.90, kind="other", anchor="end"),
    dict(name="YASPS", x=0.59, y=0.87, kind="other", anchor="start"),
    dict(name="SymX", x=0.55, y=0.64, kind="other", anchor="start"),
    dict(name="STARK", x=0.79, y=0.83, kind="other", anchor="end"),
    dict(name="PolyFEM", x=0.75, y=0.75, kind="uses", anchor="end"),
    dict(name="Rigid IPC", x=0.95, y=0.72, kind="uses", anchor="end"),
    dict(name="SfePy / FEDOO / MeshFEM", x=0.95, y=0.62, kind="uses", anchor="end"),
    # Integrated frameworks with other contact models
    dict(name="FEBio / Abaqus / Ansys / LS-DYNA", x=0.62, y=0.42, kind="other", anchor="start"),
    dict(name="SOFA", x=0.54, y=0.31, kind="other", anchor="start"),
    dict(name="Chrono / Drake", x=0.76, y=0.30, kind="other", anchor="start"),
    dict(name="Bullet / PhysX / MuJoCo", x=0.60, y=0.14, kind="other", anchor="start"),
    dict(name="NVIDIA Newton", x=0.90, y=0.20, kind="other", anchor="end"),
]

# Canvas geometry in SVG user units. The figure scales with its container.
WIDTH, HEIGHT = 780, 520
LEFT, TOP, RIGHT, BOTTOM = 64, 16, 16, 86
PLOT_W = WIDTH - LEFT - RIGHT
PLOT_H = HEIGHT - TOP - BOTTOM

RADIUS = {"toolkit": 11, "uses": 7, "other": 7}


def sx(x):
    return LEFT + x * PLOT_W


def sy(y):
    # SVG's y axis points down, so we flip it to put "IPC-grade" on top.
    return TOP + (1 - y) * PLOT_H


def build():
    out = []
    add = out.append
    add(
        f'<div class="ipctk-landscape">\n'
        f'<svg viewBox="0 0 {WIDTH} {HEIGHT}" xmlns="http://www.w3.org/2000/svg" '
        f'role="img" aria-labelledby="landscape-title landscape-desc">\n'
        f'<title id="landscape-title">Where the IPC Toolkit fits</title>\n'
        f'<desc id="landscape-desc">Projects placed by integration model '
        f"(modular component to integrated framework) and contact handling "
        f"(generic to IPC-grade). The IPC Toolkit is the only project in the "
        f"modular, IPC-grade corner.</desc>"
    )

    # Quadrants. The top-left one is the open territory we care about.
    mid_x, mid_y = sx(0.5), sy(0.5)
    add(f'<rect class="ls-quad" x="{LEFT}" y="{TOP}" width="{PLOT_W}" height="{PLOT_H}"/>')
    add(
        f'<rect class="ls-open" x="{LEFT}" y="{TOP}" '
        f'width="{mid_x - LEFT}" height="{mid_y - TOP}"/>'
    )
    add(f'<line class="ls-grid" x1="{mid_x}" y1="{TOP}" x2="{mid_x}" y2="{TOP + PLOT_H}"/>')
    add(f'<line class="ls-grid" x1="{LEFT}" y1="{mid_y}" x2="{LEFT + PLOT_W}" y2="{mid_y}"/>')
    add(
        f'<text class="ls-note" x="{LEFT + 10}" y="{TOP + 18}">'
        f"Open territory: IPC-grade contact you can embed</text>"
    )

    # Axes.
    base = TOP + PLOT_H
    add(f'<text class="ls-axis" x="{LEFT}" y="{base + 20}">MODULAR COMPONENTS</text>')
    add(
        f'<text class="ls-axis" x="{LEFT + PLOT_W}" y="{base + 20}" '
        f'text-anchor="end">INTEGRATED FRAMEWORKS</text>'
    )
    add(
        f'<text class="ls-axis-title" x="{LEFT + PLOT_W / 2}" y="{base + 42}" '
        f'text-anchor="middle">Integration model →</text>'
    )
    add(
        f'<text class="ls-axis" x="{LEFT - 10}" y="{TOP + 12}" '
        f'text-anchor="end" transform="rotate(-90 {LEFT - 10} {TOP + 12})">'
        f"IPC-GRADE</text>"
    )
    add(
        f'<text class="ls-axis" x="{LEFT - 10}" y="{base}" text-anchor="start" '
        f'transform="rotate(-90 {LEFT - 10} {base})">GENERIC</text>'
    )
    add(
        f'<text class="ls-axis-title" x="{LEFT - 40}" y="{TOP + PLOT_H / 2}" '
        f'text-anchor="middle" transform="rotate(-90 {LEFT - 40} {TOP + PLOT_H / 2})">'
        f"Contact handling: barrier + CCD →</text>"
    )

    # Points. We draw labels after dots so a label is never hidden by a dot.
    for p in POINTS:
        r = RADIUS[p["kind"]]
        add(f'<circle class="ls-dot ls-{p["kind"]}" cx="{sx(p["x"]):.1f}" cy="{sy(p["y"]):.1f}" r="{r}"/>')
    for p in POINTS:
        r = RADIUS[p["kind"]]
        gap = r + 6
        x = sx(p["x"]) + (gap if p["anchor"] == "start" else -gap)
        y = sy(p["y"]) + 4 + p.get("dy", 0)
        cls = "ls-label ls-label-toolkit" if p["kind"] == "toolkit" else "ls-label"
        add(
            f'<text class="{cls}" x="{x:.1f}" y="{y:.1f}" '
            f'text-anchor="{p["anchor"]}">{escape(p["name"])}</text>'
        )

    # Legend, tucked under the x axis title.
    ly = HEIGHT - 14
    items = [
        ("toolkit", "IPC Toolkit", LEFT + 150),
        ("uses", "Built on the IPC Toolkit", LEFT + 290),
        ("other", "Other projects", LEFT + 500),
    ]
    for kind, text, lx in items:
        add(f'<circle class="ls-dot ls-{kind}" cx="{lx}" cy="{ly - 4}" r="6"/>')
        add(f'<text class="ls-legend" x="{lx + 12}" y="{ly}">{escape(text)}</text>')

    add("</svg>\n</div>\n")
    return "\n".join(out)


def write():
    """Write ``_static/landscape.html``, but only if its contents changed.

    ``conf.py`` calls this at the start of every docs build. We skip the write
    when nothing changed so an unchanged figure does not touch the file's
    modification time and force Sphinx to re-read the page that embeds it.
    """
    target = pathlib.Path(__file__).parents[2] / "landscape.html"
    html = build()
    if not target.exists() or target.read_text() != html:
        target.write_text(html)
        return target, True
    return target, False


if __name__ == "__main__":
    target, changed = write()
    print(f"{'Wrote' if changed else 'Up to date:'} {target}")
