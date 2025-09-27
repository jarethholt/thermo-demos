from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colormaps

from optimize import find_equilibrium_state
from thermo import (
    CRITICAL_CHEM_POT,
    CRITICAL_TEMPERATURE,
    max_gas_chem_pot,
    min_liq_chem_pot,
)

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.colorbar import Colorbar
    from matplotlib.colors import Colormap
    from matplotlib.figure import Figure
    from matplotlib.image import AxesImage
    from matplotlib.lines import Line2D

MIN_TEMPERATURE = 0.5
MAX_TEMPERATURE = 1.5
NUM_TEMPERATURE_PTS = 100
MIN_CHEM_POT = -2.5
MAX_CHEM_POT = -1.5
NUM_CHEM_POT_PTS = 100
FIGSIZE = (5.5, 4.5)
COLORMAP = colormaps["Blues"]

LABEL_TEXT_KWARGS = {"fontsize": 14}
CRITICAL_LINE_KWARGS = {"color": "gray"}
METASTABLE_BOUNDARY_KWARGS = {"color": "gray", "linestyle": "--"}

Array1D = np.ndarray[tuple[int], np.dtype[np.float64]]
Array2D = np.ndarray[tuple[int, int], np.dtype[np.float64]]


def _COLORMAP_DEFAULT_FACTORY() -> "Colormap":
    return COLORMAP


@dataclass(kw_only=True)
class PhaseArray:
    min_temperature: float = MIN_TEMPERATURE
    max_temperature: float = MAX_TEMPERATURE
    num_temperature_pts: int = NUM_TEMPERATURE_PTS
    min_chem_pot: float = MIN_CHEM_POT
    max_chem_pot: float = MAX_CHEM_POT
    num_chem_pot_pts: int = NUM_CHEM_POT_PTS

    def __post_init__(self):
        self.temperature = PhaseArray._get_midpoint_array(
            self.min_temperature, self.max_temperature, self.num_temperature_pts
        )
        self.chem_pot = PhaseArray._get_midpoint_array(
            self.min_chem_pot, self.max_chem_pot, self.num_chem_pot_pts
        )
        self.density: Array2D = np.array(
            [
                [find_equilibrium_state(t, c).density for c in self.chem_pot]
                for t in self.temperature
            ]
        )

    @staticmethod
    def _get_midpoint_array(
        min_value: float, max_value: float, num_values: int
    ) -> Array1D:
        return min_value + (max_value - min_value) / num_values * (
            0.5 + np.arange(num_values, dtype=np.float64)
        )

    @staticmethod
    def _get_edge_array(min_value: float, max_value: float, num_values: int) -> Array1D:
        return min_value + (max_value - min_value) / num_values * np.arange(
            num_values + 1, dtype=np.float64
        )

    def set_up_axes(self, axes: "Axes") -> None:
        axes.set_xlim(self.min_temperature, self.max_temperature)
        axes.set_xticks(
            [self.min_temperature, CRITICAL_TEMPERATURE, self.max_temperature]
        )
        axes.set_xlabel("Temperature", **LABEL_TEXT_KWARGS)  # type: ignore
        axes.set_ylim(self.min_chem_pot, self.max_chem_pot)
        axes.set_yticks([self.min_chem_pot, CRITICAL_CHEM_POT, self.max_chem_pot])
        axes.set_ylabel("Chem pot", **LABEL_TEXT_KWARGS)  # type: ignore
        axes.set_aspect("equal")
        axes.set_title("Equilibrium density", **LABEL_TEXT_KWARGS)  # type: ignore

    def add_critical_lines(self, axes: "Axes") -> tuple["Line2D", "Line2D"]:
        (critical_chem_pot_line,) = axes.plot(
            [self.min_temperature, CRITICAL_TEMPERATURE],
            [CRITICAL_CHEM_POT, CRITICAL_CHEM_POT],
            **CRITICAL_LINE_KWARGS,  # type: ignore
        )
        (critical_temperature_line,) = axes.plot(
            [CRITICAL_TEMPERATURE, CRITICAL_TEMPERATURE],
            [self.min_chem_pot, self.max_chem_pot],
            **CRITICAL_LINE_KWARGS,  # type: ignore
        )
        return critical_chem_pot_line, critical_temperature_line

    def add_metastable_lines(self, axes: "Axes") -> tuple["Line2D", "Line2D"]:
        subcrit_temperature = PhaseArray._get_edge_array(
            self.min_temperature, CRITICAL_TEMPERATURE, self.num_chem_pot_pts // 2
        )
        (gas_metastable_boundary,) = plt.plot(
            subcrit_temperature,
            max_gas_chem_pot(subcrit_temperature),
            **METASTABLE_BOUNDARY_KWARGS,  # type: ignore
        )
        (liq_metastable_boundary,) = plt.plot(
            subcrit_temperature,
            min_liq_chem_pot(subcrit_temperature),
            **METASTABLE_BOUNDARY_KWARGS,  # type: ignore
        )
        return gas_metastable_boundary, liq_metastable_boundary

    def draw_phase_diagram(
        self, figure: "Figure", axes: "Axes", colormap: "Colormap | str"
    ) -> tuple["AxesImage", "Colorbar"]:
        image = axes.imshow(
            self.density.T,
            cmap=colormap,
            vmin=0.0,
            vmax=1.0,
            origin="lower",
            extent=(
                self.min_temperature,
                self.max_temperature,
                self.min_chem_pot,
                self.max_chem_pot,
            ),
        )
        colorbar = figure.colorbar(image, ax=axes)
        return image, colorbar


@dataclass(kw_only=True)
class PhaseDiagram:
    phase_array: PhaseArray
    figsize: tuple[float, float] = FIGSIZE
    colormap: "Colormap | str" = field(default_factory=_COLORMAP_DEFAULT_FACTORY)

    def __post_init__(self):
        self.figure, self.axes = self._initialize_figure()
        (
            self.critical_chem_pot_line,
            self.critical_temperature_line,
            self.gas_metastable_boundary,
            self.liq_metastable_boundary,
        ) = self._add_guide_lines(self.axes)
        self.image, self.colorbar = self.phase_array.draw_phase_diagram(
            self.figure, self.axes, self.colormap
        )

    def _initialize_figure(self) -> tuple["Figure", "Axes"]:
        figure, axes = plt.subplots(figsize=self.figsize)
        self.phase_array.set_up_axes(axes)
        return figure, axes

    def _add_guide_lines(
        self, axes: "Axes"
    ) -> tuple["Line2D", "Line2D", "Line2D", "Line2D"]:
        critical_chem_pot_line, critical_temperature_line = (
            self.phase_array.add_critical_lines(axes)
        )
        gas_metastable_boundary, liq_metastable_boundary = (
            self.phase_array.add_metastable_lines(axes)
        )
        return (
            critical_chem_pot_line,
            critical_temperature_line,
            gas_metastable_boundary,
            liq_metastable_boundary,
        )
