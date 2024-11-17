import numpy as np
from manim import Mobject, Animation, linear
from numpy.typing import NDArray


class MoveAlongXYZPath(Animation):
    """
    An animation that moves an Mobject along a path defined by an array of points.

    If no timestamp array is given, the path will be taken in linear time.
    Otherwise the timestamp array will be used.

    The `rate_func` parameter can still be used and will act as function composition with the timestamp array.

    When finding the position at time t, the closest value in the position array
    will be used. No interpolation is performed.
    ----------
    mobject:
        The object to be animated
    *path:
        An array of (x,y,z) points, optionally preceeded by an equal length array of timestamps
    run_time:
        The run time of the animation. (ignored if a timestamp array is provided)
    """

    def __init__(
        self,
        mobject: Mobject,
        *path: NDArray,
        run_time: float | int = 1,
        rate_func: callable = linear,
        suspend_mobject_updating: bool = False,
        **kwargs,
    ) -> None:
        match len(path):
            case 2:
                ts, points = path
                run_time = np.max(ts)
                assert len(ts) == len(points), "Arrays are not the same length"
            case 1:
                points = path[0]
                ts = np.linspace(0, run_time, len(points))
            case _:
                raise Exception("wrong number of arguments")

        # Sort the array so lookups can be performed faster
        ts, points = map(np.array, zip(*sorted([*zip(ts, points)])))

        # The `Animation class` wants a [0,1] time frame
        self.alphas = ts / run_time
        self.points = points

        super().__init__(
            mobject,
            suspend_mobject_updating=suspend_mobject_updating,
            run_time=run_time,
            rate_func=rate_func,
            **kwargs,
        )

    def interpolate_mobject(self, alpha: float) -> None:
        index = np.searchsorted(self.alphas, alpha)
        point = self.points[index]
        self.mobject.move_to(point)
