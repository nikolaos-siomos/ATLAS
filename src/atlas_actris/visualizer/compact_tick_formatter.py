"""Scientific y-axis formatting triggered by typographic tick-label width."""
import math
import numpy as np
from matplotlib.textpath import TextToPath
from matplotlib.ticker import ScalarFormatter


class CompactScalarFormatter(ScalarFormatter):
    """Switch to a shared power of ten when labels exceed a reference width.

    Width is measured in points using the tick font and Matplotlib's text/math
    metrics, so decimal points and minus signs contribute their actual widths.
    The reference text specifies a width, not a numeric threshold. Available
    space between panels can impose a smaller limit, reserving room for the
    axis title and padding. Precision,
    tick locations, axis limits and subplot geometry remain independent.
    """

    def __init__(self, max_label, *, check_panel_gap=True, max_plain_magnitude=None):
        super().__init__(useOffset=False, useMathText=False, useLocale=False)
        self.max_label = max_label
        self.check_panel_gap = check_panel_gap
        self.max_plain_magnitude = max_plain_magnitude
        self._text_metrics = TextToPath()
        self.set_powerlimits((0, 0))

    def _label_width(self, text, font):
        return self._text_metrics.get_text_width_height_descent(
            text, font, ismath=text.startswith("$")
        )[0]

    def _available_tick_width(self):
        """Space in points left after the y-title, pads and a 4 pt clearance."""
        if not self.check_panel_gap:
            return math.inf
        ax = self.axis.axes
        position = ax.get_position()
        neighbours = [other for other in ax.figure.axes if other is not ax
                      and other.get_position().x1 <= position.x0
                      and other.get_position().y1 > position.y0
                      and other.get_position().y0 < position.y1]
        if not neighbours:
            return math.inf
        previous_edge = max(other.get_position().x1 for other in neighbours)
        gap_points = (position.x0 - previous_edge) * ax.figure.get_figwidth() * 72.0
        tick = self.axis.get_major_ticks(1)[0]
        tick_padding = tick.get_pad() + tick.get_tick_padding()
        title = self.axis.label
        title_space = (1.2 * title.get_fontsize() + self.axis.labelpad) if title.get_text() else 0.0
        return max(0.0, gap_points - tick_padding - title_space - 4.0)

    def set_locs(self, locs):
        self.set_scientific(False)
        if self.axis is None:
            super().set_locs(locs)
            return
        lower, upper = sorted(self.axis.get_view_interval())
        visible = np.asarray([value for value in locs
                              if math.isfinite(value) and lower <= value <= upper])
        # Out-of-view locator ticks must not enlarge precision or label width.
        working_locs = visible if len(visible) else locs
        super().set_locs(working_locs)
        if len(visible) == 0:
            return
        font = self.axis.get_major_ticks(1)[0].label1.get_fontproperties()
        reference = self.fix_minus(self.max_label)
        if self.get_useMathText():
            reference = r"$\mathdefault{" + reference + "}$"
        limit = min(self._label_width(reference, font), self._available_tick_width())
        labels = [super(CompactScalarFormatter, self).__call__(value) for value in visible]
        too_large = self.max_plain_magnitude is not None and any(
            abs(value) >= self.max_plain_magnitude for value in visible
        )
        # Also catch a label that rounds up to the magnitude threshold.
        if self.max_plain_magnitude is not None:
            for label in labels:
                try:
                    too_large = too_large or abs(float(label.replace("−", "-"))) >= self.max_plain_magnitude
                except ValueError:
                    pass
        if too_large or any(self._label_width(label, font) > limit + 1e-6 for label in labels):
            self.set_scientific(True)
            super().set_locs(working_locs)
