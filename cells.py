from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np


@dataclass
class Points:
    """
    2D points with x and y coordinates.
    """

    array: np.ndarray

    @property
    def x(self) -> np.ndarray:
        return self.array[:, 0]

    @property
    def y(self) -> np.ndarray:
        return self.array[:, 1]

    def __len__(self):
        return self.array.shape[0]

    def __getitem__(self, item: Any) -> Points:
        return Points(array=self.array[item])

    @classmethod
    def generate(cls, n: int) -> Points:
        """
        Generate n points using a random uniform distribution.
        """
        return Points(array=np.random.rand(n, 2) - np.array([0.5, 0.5]))

    def draw(self, ax: Optional[plt.Axes] = None, *args, **kwargs) -> Any:
        """
        Draw the points on the given axes.
        """
        if ax is None:
            ax = plt.gca()
        return ax.scatter(self.x, self.y, *args, **kwargs)


class Cell(ABC):
    """
    Abstract class for a cell.
    """

    @abstractmethod
    def sample(self, n: int) -> Points:
        """
        Sample uniformly n points within the current cell.
        """
        pass

    @abstractmethod
    def draw(self, ax: Optional[plt.Axes] = None, *args, **kwargs) -> Any:
        """
        Draw the cell structure on the given axes.
        """
        pass


@dataclass
class Circle(Cell):
    """
    Circular cell.
    """

    radius: float = 1.0
    center: np.ndarray = np.array([0.0, 0.0])

    def __post_init__(self):
        self.center = np.asarray(self.center).reshape(1, 2)

    def sample(self, n: int) -> Points:
        def sample_in_circle(n):
            xy = np.random.rand(n, 2)
            xy = (xy - 0.5) * (2 * self.radius)
            r = np.linalg.norm(xy, axis=1)

            return xy[r <= self.radius, :]

        n2 = n << 1
        xy = sample_in_circle(n2)

        while xy.shape[0] < n:
            xy = sample_in_circle(n2)

        return Points(array=xy[:n, :] + self.center)

    def draw(self, ax: Optional[plt.Axes] = None, *args, fill=False, **kwargs) -> Any:
        if ax is None:
            ax = plt.gca()
        circle = plt.Circle(
            self.center.flatten(), self.radius, *args, fill=fill, **kwargs
        )
        return ax.add_patch(circle)


@dataclass
class Hexagon(Cell):
    """
    Hexagonal cell.
    """

    radius: float = 1.0
    center: np.ndarray = np.array([0.0, 0.0])
    rotation: float = 0.0

    def __post_init__(self):
        self.center = np.asarray(self.center).reshape(1, 2)

    @property
    def corners_centered_array(self) -> np.ndarray:
        zero = 0.0
        radi = self.radius
        s3_2 = np.sqrt(3) * 0.5 * self.radius
        ra_2 = self.radius * 0.5
        array = np.array(
            [
                [+s3_2, +ra_2],
                [+zero, +radi],
                [-s3_2, +ra_2],
                [-s3_2, -ra_2],
                [+zero, -radi],
                [+s3_2, -ra_2],
            ]
        )
        return array

    @property
    def rotation_array(self) -> np.ndarray:
        cos = np.cos(self.rotation)
        sin = np.sin(self.rotation)

        return np.array(
            [
                [+cos, -sin],
                [+sin, +cos],
            ]
        )

    @property
    def corners_centered(self) -> Points:
        array = self.corners_centered_array
        return Points(array=array)

    @property
    def corners_array(self) -> np.ndarray:
        array = self.corners_centered_array
        array = (self.rotation_array @ array.T).T
        return array + self.center

    @property
    def corners(self) -> Points:
        array = self.corners_array
        return Points(array=array)

    def sample(self, n: int) -> Points:
        def sample_in_hexagon(n):
            """
            We check that points are in hexagon using
            the 'cross product trick'.
            """
            xy = np.random.rand(n, 2)
            xy = (xy - 0.5) * (2 * self.radius)
            corners = self.corners_centered_array

            # Vectors from points (samples) to each corner
            v1 = -xy[:, :, None] + corners.T[None, :, :]
            v2 = -xy[:, :, None] + np.roll(corners, 1, axis=0).T[None, :, :]

            cross = np.cross(v1, v2, axis=1)

            # Sum is 6 (or -6) if all the signs are the same
            s = np.sum(np.sign(cross), axis=1)
            index = np.abs(s) == 6

            return xy[index, :]

        n2 = n << 1
        xy = sample_in_hexagon(n2)

        while xy.shape[0] < n:
            xy = sample_in_hexagon(n2)

        xy = (self.rotation_array @ xy.T).T

        return Points(array=xy[:n, :] + self.center)

    def draw(self, ax: Optional[plt.Axes] = None, *args, fill=False, **kwargs) -> Any:
        if ax is None:
            ax = plt.gca()
        hexagon = plt.Polygon(self.corners_array, *args, fill=fill, **kwargs)
        return ax.add_patch(hexagon)
