import numpy as np
import xarray as xr


def _to_numpy_copy(a, name):
    if isinstance(a, xr.DataArray):
        return a.copy(deep=True).values
    if isinstance(a, np.ndarray):
        return np.copy(a)
    raise TypeError(f"Unsupported type for {name}: {type(a)}")


def _check_err_type(err_type):
    if err_type not in ["sem", "std"]:
        raise ValueError(
            "Error type provided is not understood. "
            "Please select one of: sem, std"
        )


def _as_odd_window(win):
    win = int(win)

    if win < 1:
        win = 1

    if win % 2 == 0:
        win += 1

    return win


def _coord_step(vals):
    vals = np.asarray(vals)

    if vals.size < 2:
        raise ValueError("Coordinate array must contain at least two points.")

    return vals[1] - vals[0]


def _limits_to_slice(vals, lims):
    """
    Convert physical limits to a clipped half-open slice [s_bin:e_bin).

    Limits outside the coordinate range are allowed.
    If there is no overlap with the signal, an empty slice is returned.
    """

    vals = np.asarray(vals)

    lo = min(lims[0], lims[-1])
    hi = max(lims[0], lims[-1])

    mask = (vals >= lo) & (vals <= hi)

    if not np.any(mask):
        return 0, 0

    inds = np.where(mask)[0]

    return inds[0], inds[-1] + 1

def _rolling_mean_std_1d_left_partial_right_required(y, win):
    y = np.asarray(y, dtype=float)
    win = _as_odd_window(win)
    buf = win // 2

    avg = np.full(y.shape, np.nan, dtype=float)
    std = np.full(y.shape, np.nan, dtype=float)
    count = np.zeros(y.shape, dtype=float)

    n_bins = y.size

    finite = np.isfinite(y)
    y0 = np.where(finite, y, 0.0)

    c_count = np.concatenate(([0.0], np.cumsum(finite.astype(float))))
    c_sum = np.concatenate(([0.0], np.cumsum(y0)))
    c_sum2 = np.concatenate(([0.0], np.cumsum(y0 ** 2)))

    for i in range(n_bins):
        lo = max(0, i - buf)
        hi = i + buf + 1

        if hi > n_bins:
            continue

        n = c_count[hi] - c_count[lo]

        count[i] = n

        if n == 0:
            continue

        s = c_sum[hi] - c_sum[lo]
        s2 = c_sum2[hi] - c_sum2[lo]

        mean = s / n
        var = s2 / n - mean ** 2

        if var < 0.0:
            var = 0.0

        avg[i] = mean
        std[i] = np.sqrt(var)

    return avg, std, count


def _rolling_mean_std_2d_left_partial_right_required(z, win):
    z = np.asarray(z, dtype=float)
    win = _as_odd_window(win)
    buf = win // 2

    avg = np.full(z.shape, np.nan, dtype=float)
    std = np.full(z.shape, np.nan, dtype=float)
    count = np.zeros(z.shape, dtype=float)

    n_bins = z.shape[1]

    finite = np.isfinite(z)
    z0 = np.where(finite, z, 0.0)

    c_count = np.concatenate(
        [np.zeros((z.shape[0], 1)), np.cumsum(finite.astype(float), axis=1)],
        axis=1,
    )
    c_sum = np.concatenate(
        [np.zeros((z.shape[0], 1)), np.cumsum(z0, axis=1)],
        axis=1,
    )
    c_sum2 = np.concatenate(
        [np.zeros((z.shape[0], 1)), np.cumsum(z0 ** 2, axis=1)],
        axis=1,
    )

    for i in range(n_bins):
        lo = max(0, i - buf)
        hi = i + buf + 1

        if hi > n_bins:
            continue

        n = c_count[:, hi] - c_count[:, lo]

        count[:, i] = n

        valid = n > 0

        if not np.any(valid):
            continue

        s = c_sum[:, hi] - c_sum[:, lo]
        s2 = c_sum2[:, hi] - c_sum2[:, lo]

        mean = np.full(z.shape[0], np.nan, dtype=float)
        var = np.full(z.shape[0], np.nan, dtype=float)

        mean[valid] = s[valid] / n[valid]
        var[valid] = s2[valid] / n[valid] - mean[valid] ** 2
        var[valid] = np.maximum(var[valid], 0.0)

        avg[:, i] = mean
        std[:, i] = np.sqrt(var)

    return avg, std, count


def sliding_average_1D_fast(
        y_vals, x_vals, x_sm_lims, x_sm_win,
        expo=None, err_type="sem"):

    _check_err_type(err_type)

    x_vals = _to_numpy_copy(x_vals, "x_vals")
    y_vals = _to_numpy_copy(y_vals, "y_vals").astype(float)

    dx = _coord_step(x_vals)

    win = _as_odd_window(1E-3 * x_sm_win / dx)

    s_bin, e_bin = _limits_to_slice(x_vals, x_sm_lims)

    y_avg = np.copy(y_vals).astype(float)
    y_err = np.full(y_vals.shape, np.nan, dtype=float)

    if s_bin >= e_bin:
        return y_avg, y_err

    avg_all, std_all, count_all = _rolling_mean_std_1d_left_partial_right_required(
        y_vals,
        win=win,
    )

    y_avg[s_bin:e_bin] = avg_all[s_bin:e_bin]

    if err_type == "sem":
        with np.errstate(divide="ignore", invalid="ignore"):
            y_err[s_bin:e_bin] = (
                std_all[s_bin:e_bin] / np.sqrt(count_all[s_bin:e_bin])
            )
    elif err_type == "std":
        y_err[s_bin:e_bin] = std_all[s_bin:e_bin]

    return y_avg, y_err


def sliding_average_2D_fast(
        z_vals, y_vals, y_sm_lims, y_sm_win,
        expo=None, err_type="sem"):

    _check_err_type(err_type)

    z_vals = _to_numpy_copy(z_vals, "z_vals").astype(float)
    y_vals = _to_numpy_copy(y_vals, "y_vals")

    dy = _coord_step(y_vals)

    win = _as_odd_window(1E-3 * y_sm_win / dy)

    s_bin, e_bin = _limits_to_slice(y_vals, y_sm_lims)

    z_avg = np.copy(z_vals).astype(float)
    z_err = np.full(z_vals.shape, np.nan, dtype=float)

    if s_bin >= e_bin:
        return z_avg, z_err

    avg_all, std_all, count_all = _rolling_mean_std_2d_left_partial_right_required(
        z_vals,
        win=win,
    )

    z_avg[:, s_bin:e_bin] = avg_all[:, s_bin:e_bin]

    if err_type == "sem":
        with np.errstate(divide="ignore", invalid="ignore"):
            z_err[:, s_bin:e_bin] = (
                std_all[:, s_bin:e_bin] / np.sqrt(count_all[:, s_bin:e_bin])
            )
    elif err_type == "std":
        z_err[:, s_bin:e_bin] = std_all[:, s_bin:e_bin]

    return z_avg, z_err


def sliding_average_2D_bin_fast(
        z_vals, y_vals, y_sm_lims, y_sm_win,
        expo=None, err_type="sem"):

    _check_err_type(err_type)

    z_vals = _to_numpy_copy(z_vals, "z_vals").astype(float)
    y_vals = _to_numpy_copy(y_vals, "y_vals")

    win = _as_odd_window(y_sm_win)

    s_bin = int(y_sm_lims[0])
    e_bin = int(y_sm_lims[1])

    s_bin = max(s_bin, 0)
    e_bin = min(e_bin, z_vals.shape[1])

    z_avg = np.copy(z_vals).astype(float)
    z_err = np.full(z_vals.shape, np.nan, dtype=float)

    if s_bin >= e_bin:
        return z_avg, z_err

    avg_all, std_all, count_all = _rolling_mean_std_2d_left_partial_right_required(
        z_vals,
        win=win,
    )

    z_avg[:, s_bin:e_bin] = avg_all[:, s_bin:e_bin]

    if err_type == "sem":
        with np.errstate(divide="ignore", invalid="ignore"):
            z_err[:, s_bin:e_bin] = (
                std_all[:, s_bin:e_bin] / np.sqrt(count_all[:, s_bin:e_bin])
            )
    elif err_type == "std":
        z_err[:, s_bin:e_bin] = std_all[:, s_bin:e_bin]

    return z_avg, z_err