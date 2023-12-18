from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np
from cells import Points


class Obstacles(ABC):
    """
    Abstract class for a list of obstacles.
    """

    @abstractmethod
    def contains(self, points: Points) -> np.ndarray:
        """
        Return a boolean array with the same size
        as the number of point.

        An entry is True if the corresponding point is
        within as least one of the obstacles.
        """
        pass

    @abstractmethod
    def draw(self, ax: Optional[plt.Axes] = None, *args, **kwargs):
        """
        Draw the obstacles on the given axes.
        """
        pass


@dataclass
class CircularObstacles(Obstacles):
    """
    Circular cell.
    """

    radii: float = np.array([0.1])
    centers: np.ndarray = np.array([0.0, 0.0])

    def __post_init__(self):
        self.radii = np.asarray(self.radii).reshape(-1)
        self.centers = np.asarray(self.centers).reshape(-1, 2)

    @classmethod
    def generate(
        cls, n: int, radius: float = 0.1, radius_std: float = 0.0
    ) -> CircularObstacles:
        """
        Generate n, possibly-overlapping, obstacles,
        using a random uniform distribution for
        the centers, and a Gaussian distribution for the radii.
        """
        rng = np.random.default_rng()
        radii = rng.normal(loc=radius, scale=radius_std, size=n)
        centers = np.random.rand(n, 2) - np.array([0.5, 0.5])

        return CircularObstacles(radii=radii, centers=centers)

    def contains(self, points: np.ndarray) -> np.ndarray:
        xp, yp = points.T
        xc, yc = self.centers.T
        dx = np.subtract.outer(xp, xc)
        dy = np.subtract.outer(yp, yc)
        d2 = dx * dx + dy * dy
        r2 = self.radii * self.radii

        cond = d2 <= r2

        return np.any(cond, axis=1)

    def draw(self, ax: Optional[plt.Axes] = None, *args, fill=False, **kwargs):
        if ax is None:
            ax = plt.gca()

        artists = []
        for radius, center in zip(self.radii, self.centers):
            circle = plt.Circle(center.flatten(), radius, *args, fill=fill, **kwargs)
            artists.append(ax.add_patch(circle))

        return artists
