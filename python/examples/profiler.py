import pandas as pd
import plotly.express as px
import pathlib

def plot(csv_path):
    df = pd.read_csv(csv_path, na_filter=False)
    df["Time (s)"] = df["Time (ms)"] / 1000.0
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

plot(pathlib.Path(__file__).parent / "lbvh_profile.csv")