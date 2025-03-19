import matplotlib.pyplot as plt
import numpy as np
import argparse


def get_ref_sol(opts):
    ref_sol = np.loadtxt(f"data/ref_{opts.grid_pts_1d}.txt")
    print(f"Reference solution generated in {ref_sol[0]}s")

    if opts.plot:
        u = ref_sol[1::2].reshape(opts.grid_pts_1d, opts.grid_pts_1d)
        v = ref_sol[2::2].reshape(opts.grid_pts_1d, opts.grid_pts_1d)

        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        extent = [-1, 1, -1, 1]
        axes[0].imshow(u, cmap="viridis", extent=extent)
        axes[0].set_title("u")

        axes[1].imshow(v, cmap="viridis", extent=extent)
        axes[1].set_title("v")

        plt.show()

    return ref_sol


def sol_error(ref_sol, sol):
    return np.linalg.norm(ref_sol[1:] - sol[1:]) / np.linalg.norm(ref_sol[1:])


def get_lsrk_results(ref_sol, opts):
    method_names = ("LSRK", "ERK", "DIRK")
    tols = (f"1e-{i}" for i in range(2, 10))
    times = np.zeros((len(tols), len(method_names)))
    errors = np.zeros_like(times)

    for tol_index in range(len(tols)):
        tol = tols[tol_index]
        for method_index in range(len(method_names)):
            method_name = method_names[method_index]
            sol = np.loadtxt(f"data/{method_name}_{tol}.txt")
            times[tol_index, method_index] = sol[0]
            errors[tol_index, method_index] = sol_error(ref_sol, sol)

    header = (f"{n} {t}" for n in method_names for t in ("Time", "Error"))
    np.savetxt(
        "data/lsrk.csv", np.hstack((times, errors)), header=", ".join(header)
    )

    if opts.plot:
        plt.loglog(times, errors)
        plt.title("LSRK Performance")
        plt.xlabel("Time (s)")
        plt.ylabel("Relative l2 Error")
        plt.show()


def get_splitting_results(ref_sol, opts):
    orders = (1, 2, 3, 4, 6)
    dts = 2 ** np.arange(-10, -4)
    errors = np.zeros((len(dts), len(orders)))

    for dt_index in range(len(dts)):
        dt = dts[dt_index]
        for order_index in range(len(orders)):
            order = orders[order_index]
            sol = np.loadtxt(f"data/Splitting_{order}_{dt}.txt")
            errors[dt_index, order_index] = sol_error(ref_sol, sol)

    headers = ("dt",) + (f"Order {o} Error" for o in orders)
    np.savetxt(
        "data/splitting.csv",
        np.hstack((dts, errors)),
        header=", ".join(headers),
    )

    if opts.plot:
        plt.loglog(dts, errors)
        plt.title("Splitting Method Convergence")
        plt.xlabel("Time Step")
        plt.ylabel("Relative l2 Error")
        plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot", action="store_true", default=False)
    parser.add_argument("--grid_pts_1d", type=int, default=1024)
    opts = parser.parse_args()

    ref_sol = get_ref_sol(opts)
    get_lsrk_results(ref_sol, opts)


if __name__ == "__main__":
    main()
