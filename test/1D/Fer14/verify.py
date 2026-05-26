import argparse
import sys
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

warnings.filterwarnings("ignore")

# Configure matplotlib for transparent SVG with theme-aware styling
mpl.rcParams.update(
    {
        "figure.facecolor": "none",
        "axes.facecolor": "none",
        "savefig.facecolor": "none",
        "svg.fonttype": "none",
        "font.size": 14,
        "axes.titlesize": 16,
        "axes.labelsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 12,
        "figure.titlesize": 16,
    }
)


def _setup_orion_import() -> None:
    """Add local ORION Python source to sys.path if available."""
    root = Path(__file__).resolve().parent
    for parent in [root, *root.parents]:
        candidate = parent / "lib" / "ORION" / "src" / "python"
        if candidate.exists():
            sys.path.insert(0, str(candidate))
            return


_setup_orion_import()

try:
    from ORION import read_TEC
except ModuleNotFoundError as exc:
    raise ModuleNotFoundError(
        "ORION package not found. Activate the correct environment or install ORION."
    ) from exc


def _cell_centers(x_nodes: np.ndarray) -> np.ndarray:
    """Return 1D cell-center coordinates from x nodes along the i-direction."""
    x_nodes_arr = np.asarray(x_nodes)

    if x_nodes_arr.ndim == 3:
        x_nodes_1d = x_nodes_arr[:, 0, 0]
    elif x_nodes_arr.ndim == 2:
        x_nodes_1d = x_nodes_arr[:, 0]
    else:
        x_nodes_1d = x_nodes_arr.reshape(-1)

    return 0.5 * (x_nodes_1d[:-1] + x_nodes_1d[1:])


def _extract_profiles(file_path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read one TEC file and extract x cell centers, x-velocity, and temperature."""
    x_arr, _y_arr, _z_arr, var_arr, _vnames = read_TEC(str(file_path))

    x_centers = _cell_centers(x_arr[0])
    velocity_x = np.asarray(var_arr[0][9]).reshape(-1)
    temperature = np.asarray(var_arr[0][13]).reshape(-1)

    n = min(x_centers.size, velocity_x.size, temperature.size)
    return x_centers[:n], velocity_x[:n], temperature[:n]


def _load_reference_series(csv_path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Load reference velocity points grouped by series label."""
    ref_raw = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")

    if ref_raw.size == 0:
        return {}
    if ref_raw.shape == ():
        ref_raw = np.array([ref_raw], dtype=ref_raw.dtype)

    ref_data: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for series_name in np.unique(ref_raw["series"]):
        mask = ref_raw["series"] == series_name
        x_ref = np.asarray(ref_raw["x_m"][mask], dtype=float)
        u_ref = np.asarray(ref_raw["velocity_m_s"][mask], dtype=float)
        ref_data[str(series_name).strip().lower()] = (x_ref, u_ref)

    return ref_data


def _load_reference_series_with_y(
    csv_path: Path, y_candidates: list[str]
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Load reference points grouped by series label, selecting first available y column."""
    ref_raw = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None, encoding="utf-8")

    if ref_raw.size == 0:
        return {}
    if ref_raw.shape == ():
        ref_raw = np.array([ref_raw], dtype=ref_raw.dtype)

    names_set = set(ref_raw.dtype.names or ())
    y_col = next((c for c in y_candidates if c in names_set), None)
    if y_col is None:
        raise ValueError(
            f"Could not find any y-column in {csv_path}. Expected one of: {y_candidates}"
        )

    ref_data: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for series_name in np.unique(ref_raw["series"]):
        mask = ref_raw["series"] == series_name
        x_ref = np.asarray(ref_raw["x_m"][mask], dtype=float)
        y_ref = np.asarray(ref_raw[y_col][mask], dtype=float)
        ref_data[str(series_name).strip().lower()] = (x_ref, y_ref)

    return ref_data


def _map_fields_to_reference_keys(
    field_ids: list[int], ref_series: dict[str, tuple[np.ndarray, np.ndarray]]
) -> dict[int, str]:
    """Return field->reference mapping for both legacy and numeric labels."""
    keys = set(ref_series.keys())

    # Legacy CSVs used color names.
    color_map = {17: "blue", 19: "cyan", 23: "orange"}
    if set(color_map.values()).issubset(keys):
        return color_map

    # New CSVs use numeric labels (1,2,3): map by sorted series id.
    numeric_keys = sorted([k for k in keys if k.isdigit()], key=int)
    if len(numeric_keys) >= len(field_ids):
        return {fid: numeric_keys[i] for i, fid in enumerate(field_ids)}

    return {}


def _legend_handles() -> list[Line2D]:
    """Build a compact legend: style semantics + color-time mapping."""
    return [
        Line2D([0], [0], color="black", lw=2.5, linestyle="-", label="CFD"),
        Line2D(
            [0],
            [0],
            color="black",
            marker="o",
            linestyle="None",
            markersize=8,
            markerfacecolor="black",
            label="Reference",
        ),
        Line2D([0], [0], color="red", lw=2.5, linestyle="-", label="170 $\\mu$s"),
        Line2D([0], [0], color="green", lw=2.5, linestyle="-", label="190 $\\mu$s"),
        Line2D([0], [0], color="blue", lw=2.5, linestyle="-", label="230 $\\mu$s"),
    ]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot temperature and velocity from OUTPUT/field17.tec, field19.tec, field23.tec"
    )
    parser.add_argument(
        "--output-dir",
        default="OUTPUT",
        help="Folder containing field*.tec files (default: OUTPUT)",
    )
    parser.add_argument(
        "--save",
        default="Fer14.svg",
        help="Output plot path (default: Fer14.svg)",
    )
    parser.add_argument(
        "--reference-csv",
        default="reference/velocity.csv",
        help="Reference CSV file to overlay (default: reference/velocity.csv)",
    )
    parser.add_argument(
        "--reference-temp-csv",
        default="reference/temperature.csv",
        help="Temperature reference CSV file to overlay (default: reference/temperature.csv)",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Show interactive plot window in addition to saving the figure",
    )
    parser.add_argument(
        "--x-unit",
        choices=["cm", "m"],
        default="cm",
        help="Unit for x-axis in the plot (default: cm)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    field_ids = [17, 19, 23]

    data: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for fid in field_ids:
        field_path = output_dir / f"field{fid}.tec"
        if not field_path.exists():
            raise FileNotFoundError(f"Missing file: {field_path}")
        data[fid] = _extract_profiles(field_path)

    reference_csv = Path(args.reference_csv)
    ref_series = _load_reference_series(reference_csv) if reference_csv.exists() else {}
    field_to_ref_series = _map_fields_to_reference_keys(field_ids, ref_series)

    reference_temp_csv = Path(args.reference_temp_csv)
    ref_temp_series = (
        _load_reference_series_with_y(reference_temp_csv, ["temperature_K", "vtemperature_K"])
        if reference_temp_csv.exists()
        else {}
    )
    field_to_ref_temp_series = _map_fields_to_reference_keys(field_ids, ref_temp_series)

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))
    ax_t = axs[0]
    ax_u = axs[1]
    colors = {17: "red", 19: "green", 23: "blue"}
    x_scale = 100.0 if args.x_unit == "cm" else 1.0

    for fid in field_ids:
        x, u, t = data[fid]
        x_plot = x * x_scale
        ax_t.plot(x_plot, t, label=f"CFD field{fid}", color=colors[fid], linewidth=2)
        ax_u.plot(x_plot, u, label=f"CFD field{fid}", color=colors[fid], linewidth=2)

        ref_t_key = field_to_ref_temp_series.get(fid)
        if ref_t_key is not None and ref_t_key in ref_temp_series:
            x_t_ref, t_ref = ref_temp_series[ref_t_key]
            ax_t.scatter(
                x_t_ref * x_scale,
                t_ref,
                label=f"Reference {ref_t_key}",
                color=colors[fid],
                marker="o",
                s=20,
                facecolors=colors[fid],
                edgecolors=colors[fid],
                linewidths=5.0,
            )

        ref_key = field_to_ref_series.get(fid)
        if ref_key is not None and ref_key in ref_series:
            x_ref, u_ref = ref_series[ref_key]
            ax_u.scatter(
                x_ref * x_scale,
                u_ref,
                label=f"Reference {ref_key}",
                color=colors[fid],
                marker="o",
                s=20,
                facecolors=colors[fid],
                edgecolors=colors[fid],
                linewidths=5.0,
            )

    ax_t.set_xlabel(f"Length, {args.x_unit}", fontsize=20)
    ax_t.set_ylabel("Temperature, K", fontsize=20)
    ax_t.tick_params(axis="both", which="major", labelsize=18)
    ax_t.set_xlim(0, 12)
    ax_t.set_ylim(500, 3000)
    ax_t.grid(True, alpha=0.3)

    ax_u.set_xlabel(f"Length, {args.x_unit}", fontsize=20)
    ax_u.set_ylabel("Velocity, m/s", fontsize=20)
    ax_u.tick_params(axis="both", which="major", labelsize=18)
    ax_u.set_xlim(0, 12)
    ax_u.set_ylim(-600, 600)
    ax_u.grid(True, alpha=0.3)

    fig.legend(
        handles=_legend_handles(),
        loc="upper center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 1.03),
        columnspacing=1.2,
        handletextpad=0.6,
        fontsize=20,
    )

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])

    save_path = Path(args.save)
    fig.savefig(save_path)
    print(f"Saved figure to: {save_path}")

    if args.plot:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()