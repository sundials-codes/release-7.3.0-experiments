import matplotlib.pyplot as plt
import numpy as np

# Step Size: [
#   [||y||, ||dg/dy0||, ||dg/dp||] # order 3
#   [||y||, ||dg/dy0||, ||dg/dp||] # order 4
#   [||y||, ||dg/dy0||, ||dg/dp||] # order 5
# ]

step_sizes = [2.0**-i for i in range(1, 7)]

arkode_results = {
    3: [
        [3.4460721157298240e00, 1.7022038421161767e01, 4.2512634680695626e01],
        [1.4277626494190550e00, 3.9305067593620491e-01, 7.2034654413397281e-01],
        [1.3743209436710748e00, 2.9523474644399383e-01, 6.7990437330117570e-01],
        [1.3716318659662541e00, 2.9866295566390444e-01, 6.9908494819093681e-01],
        [1.3714766433766770e00, 2.9974295830965342e-01, 7.0256232105543714e-01],
        [1.3714674653466032e00, 2.9990809842125832e-01, 7.0304342778175255e-01],
    ],
    4: [
        [1.7225114896501261e00, 1.4612769985055976e-01, 3.2821887698773777e-01],
        [1.3787892754750568e00, 2.7637514774306371e-01, 6.3752561776181838e-01],
        [1.3716397362118002e00, 2.9696188327358297e-01, 6.9629634223171921e-01],
        [1.3714705549489457e00, 2.9971470294046060e-01, 7.0262741501122528e-01],
        [1.3714669368962780e00, 2.9991898642249959e-01, 7.0308265415413151e-01],
        [1.3714668907918908e00, 2.9993262096098472e-01, 7.0311273442238753e-01],
    ],
    5: [
        [1.4011397567944937e00, 2.0019118100388813e01, 2.3220963400067959e01],
        [1.3692127271625336e00, 3.0340002814106315e-01, 7.1445953295961062e-01],
        [1.3714480854467881e00, 3.0000883773837317e-01, 7.0331913641024535e-01],
        [1.3714666158322963e00, 2.9993474428284611e-01, 7.0311805381384984e-01],
        [1.3714668885140773e00, 2.9993356867935339e-01, 7.0311482840254780e-01],
        [1.3714668932524163e00, 2.9993355792794196e-01, 7.0311479324611770e-01],
    ],
}

julia_results = {
    3: [
        [3.4460721157298213e00, 1.7022038421162186e01, 4.2512634680696536e01],
        [1.4277626494190556e00, 3.9305067593620624e-01, 7.2034654413397436e-01],
        [1.3743209436710750e00, 2.9523474644399611e-01, 6.7990437330117948e-01],
        [1.3716318659662541e00, 2.9866295566390133e-01, 6.9908494819093003e-01],
        [1.3714766433766765e00, 2.9974295830965370e-01, 7.0256232105543859e-01],
        [1.3714674653466032e00, 2.9990809842126731e-01, 7.0304342778177120e-01],
    ]
}

# generated with OrdinaryDiffEq.jl (adaptive, tight tolerances) and Zygote.jl
reference_sol = {
    # BS3
    3: np.array([1.3714668933552083, 0.29993355811025496, 0.70311479341502030]),
    # RK4
    4: np.array([1.3714668933550507, 0.29993355811233674, 0.70311479341968513]),
    # Tsit5
    5: np.array([1.3714668933550890, 0.29993355811194400, 0.70311479341903290]),
}

# Compute absolute relative difference (entrywise)
errors = {}
for order, value_set in arkode_results.items():
    errors[order] = []
    for j, values in enumerate(value_set):
        time_step_errors = []
        for i, ref in enumerate(reference_sol[order]):
            time_step_errors.append(abs((values[i] - ref) / ref))
        errors[order].append(time_step_errors)

for order, order_error in errors.items():
    arr = np.concatenate(
        (np.reshape(step_sizes, (len(step_sizes), 1)), np.array(order_error)), 1
    )
    np.savetxt(
        f"error_asa_order_{order}.csv",
        arr,
        header="h,y,dgdy_0,dgdp",
        comments="",
        delimiter=",",
        fmt="%.16e",
    )

# Plot the error for each element
components = ["y", "dgdy_0", "dgdp"]
colors = ["r", "b", "g"]
for i, component in enumerate(components):
    plt.figure()

    j = 0
    for order, time_step_errors in errors.items():
        time_step_errors = np.array(time_step_errors)
        plt.loglog(
            step_sizes,
            time_step_errors[:, i],
            marker="o",
            color=colors[j],
            label="Error",
        )

        slope_reference = [(step / step_sizes[0]) ** order for step in step_sizes]

        plt.loglog(
            step_sizes,
            slope_reference,
            linestyle="--",
            color=colors[j],
            label=f"{order}th Order Convergence",
        )

        j = j + 1

    plt.xlabel("Time Step Size")
    plt.ylabel(f"Absolute Relative Error {component}")
    plt.title(f"Error vs Time Step Size {component}")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    plt.savefig(f"convergence_{component}.png")

# Plot the difference between arkode_results and julia_results for each order and step size
for order in arkode_results:
    if order in julia_results:
        ark = np.array(arkode_results[order])
        jul = np.array(julia_results[order])
        diff = np.abs(ark - jul)
        for i, component in enumerate(components):
            plt.figure()
            plt.plot(
                step_sizes[:len(diff)],
                diff[:, i],
                marker="o",
                color=colors[i]
            )
            plt.xlabel("Time Step Size")
            plt.ylabel(f"Abs Difference {component}")
            plt.title(
                f"ARKode vs Julia Abs Diff (order {order})\n{component}"
            )
            plt.grid(True, which="both", linestyle="--", linewidth=0.5)
            plt.savefig(
                f"arkode_vs_julia_diff_order_{order}_{component}.png"
            )
            plt.close()
