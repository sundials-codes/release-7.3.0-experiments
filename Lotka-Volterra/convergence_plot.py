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
        [3.4460721157298240e-00, 1.7022038421161767e01, 4.2512634680695626e01],
        [1.4277626494190550e-00, 3.9305067593620491e-01, 7.2034654413397281e-01],
        [1.3743209436710748e-00, 2.9523474644399383e-01, 6.7990437330117570e-01],
        [1.3716318659662541e-00, 2.9866295566390444e-01, 6.9908494819093681e-01],
        [1.3714766433766770e-00, 2.9974295830965342e-01, 7.0256232105543714e-01],
        [1.3714674653466032e-00, 2.9990809842125832e-01, 7.0304342778175255e-01],
    ],
    4: [
        [1.7061144689748420e-00, 1.0853671723496174e-00, 2.4185212896195845e-00],
        [1.3761377149646665e-00, 2.7282853937739876e-01, 6.3503646203510244e-01],
        [1.3713663679650461e-00, 2.9743750388742757e-01, 6.9789537413956060e-01],
        [1.3714544397523314e-00, 2.9976241379540819e-01, 7.0276493582762145e-01],
        [1.3714659812078640e-00, 2.9992256055665589e-01, 7.0309246082093435e-01],
        [1.3714668329282638e-00, 2.9993286237742001e-01, 7.0311338365727893e-01],
    ],
    5: [
        [1.4011397567944937e-00, 2.0019118100388813e-01, 2.3220963400067959e-01],
        [1.3692127271625336e-00, 3.0340002814106315e-01, 7.1445953295961062e-01],
        [1.3714480854467881e-00, 3.0000883773837317e-01, 7.0331913641024535e-01],
        [1.3714666158322963e-00, 2.9993474428284611e-01, 7.0311805381384984e-01],
        [1.3714668885140773e-00, 2.9993356867935339e-01, 7.0311482840254780e-01],
        [1.3714668932524163e-00, 2.9993355792794196e-01, 7.0311479324611770e-01],
    ],
}

reference_sol = np.array([1.3714668933550804, 0.2999335581120665, 0.7031147934193194])

# Compute absolute relative difference (entrywise)
errors = {}
for order, value_set in arkode_results.items():
    errors[order] = []
    for j, values in enumerate(value_set):
        time_step_errors = []
        for i, ref in enumerate(reference_sol):
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
