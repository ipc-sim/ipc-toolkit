import pandas as pd
import plotly.express as px
import pathlib
from io import StringIO

from find_ipctk import ipctk

def plot_profiler(title=None):
    df = pd.read_csv(StringIO(ipctk.profiler().csv), na_filter=False)
    # df["Time (s)"] = df["Time (ms)"] / 1000.0
    # df["Parent"] = df["Parent"].astype(str)
    # df["Id"] = df["Id"].astype(str)

    fig = px.sunburst(
        df,
        ids="Id",
        names="Name",
        parents="Parent",
        values="Time (ms)",
        branchvalues="total",
        color="Parent",
    )
    fig.update_traces(
        hovertemplate="<b>%{label}</b><br>Time (ms): %{value}<br>Call Count: %{customdata[0]}<br>Percentage of Parent: %{percentParent:.2%}",
        customdata=df[["Count"]].values,
        # tiling=dict(
        #     orientation='v'
        # )
    )
    fig.update_layout(
        title=title or "Profiler Results",
        title_x=0.5,
        title_y=0.95,
        margin=dict(t=0, l=0, r=0, b=0),
        width=800,
        height=800,
        template="plotly_white",
    )
    # fig.write_image(f"icicle_{pathlib.Path(csv_path).stem}.png", scale=2)
    # fig.write_image(f"sunburst_{pathlib.Path(csv_path).stem}.png", scale=2)
    fig.show()
    return fig

if __name__ == "__main__":
    # plot(pathlib.Path(__file__).parent / "lbvh_profile.csv")

    import meshio
    import argparse
    import numpy as np

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mesh",
        type=pathlib.Path,
        default=(pathlib.Path(__file__).parents[2] / "tests/data/puffer-ball/20.ply"))
    args = parser.parse_args()

    mesh = meshio.read(args.mesh)
    faces = mesh.cells_dict["triangle"]

    # edges = ipctk.edges(faces)
    # indices = np.lexsort((edges[:, 1], edges[:, 0]))
    # edges = edges[indices]

    # indices = np.lexsort((faces[:, 2], faces[:, 1], faces[:, 0]))
    # faces = faces[indices]

    lbvh = ipctk.LBVH()
    lbvh.build(mesh.points, edges, faces)

    plot_profiler(title=args.mesh.parent.name + "/" + args.mesh.name)
