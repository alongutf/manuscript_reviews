import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

class PanelFigure:
    """
    Helper for assembling multi‑panel figures with precise, normalized
    positioning and automatic panel labels (A, B, C…).
    """
    def __init__(self, figsize=(8, 6), label_offset=(0, 0.02),
                 label_style=None):
        self.fig = plt.figure(figsize=figsize)
        self.next_label = ord("A")
        self.label_offset = label_offset
        default_style = dict(fontweight="bold", fontsize=12,
                             ha="left", va="top")
        self.label_style = default_style or (label_style or {})

    def _auto_label(self):
        lbl = chr(self.next_label)
        self.next_label += 1
        return lbl

    def add_panel(self, rect, label=None, draw_func=None,
                  image=None, hide_axis=False):
        """
        Adds a single‑Axes panel.

        Parameters
        ----------
        rect : [l, b, w, h] in figure normalized coords.
        label : str | None
        draw_func : callable(ax) | None
        image : array | PIL image | None
        hide_axis : bool
        """
        ax = self.fig.add_axes(rect)
        if hide_axis:
            ax.set_axis_off()

        if draw_func is not None:
            draw_func(ax)

        if image is not None:
            ax.imshow(image)
            ax.set_axis_off()

        self._label_at(rect, label)
        return ax

    def add_grid_panel(self, rect, nrows, ncols,
                       label=None, sharex=False, sharey=False,
                       wspace=0.3, hspace=0.3, **subplot_kw):
        """
        Adds a panel that itself contains an nrows×ncols grid of sub‑plots.
        Returns the array of axes for further customisation.

        Everything remains vector when saved as PDF/SVG.

        Parameters
        ----------
        rect : [l, b, w, h] in figure normalized coords.
        sharex, sharey, wspace, hspace : as in plt.subplots.
        subplot_kw : forwarded to fig.add_subplot
        """
        # Convert rect (fig coords) → absolute [left, right, bottom, top]
        left, bottom, width, height = rect
        right = left + width
        top = bottom + height

        gs = GridSpec(nrows, ncols, left=left, right=right,
                      bottom=bottom, top=top,
                      wspace=wspace, hspace=hspace, figure=self.fig)

        axes = np.empty((nrows, ncols), dtype=object)
        for r in range(nrows):
            for c in range(ncols):
                axes[r, c] = self.fig.add_subplot(gs[r, c], **subplot_kw)

        # Use the first axes for placing the panel label
        self._label_at(rect, label)
        return axes

    # ------------------------------------------------------------------
    def _label_axes(self, ax, label):
        lbl = label or self._auto_label()
        ax.text(self.label_offset[0],
                1 + self.label_offset[1],
                lbl, transform=ax.transAxes,
                **self.label_style)

    def _label_at(self, rect, label):
        """
        Places the panel label at the top-left corner of the panel rectangle,
        using figure coordinates.

        Parameters
        ----------
        rect : list or tuple
            [left, bottom, width, height] of the panel in figure coordinates.
        label : str or None
            Label text (A, B, etc.). If None, auto-labels.
        """
        lbl = label or self._auto_label()
        x = rect[0] + self.label_offset[0]  # left edge + offset
        y = rect[1] + rect[3] + self.label_offset[1]  # top edge + offset
        self.fig.text(x, y, lbl, transform=self.fig.transFigure,
                      **self.label_style)

    def save(self, filename, **kwargs):
        self.fig.canvas.draw()  # ensure renderer is created
        self.fig.savefig(filename, **kwargs)