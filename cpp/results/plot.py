import csvp
import matplotlib.pyplot as plt


hdr, x, y = csvp.parse("benchmark.out"," ",float)

y *= 1e6

x = (x ** 3)[:,0]

n = 80 * x

plt.scatter(
    n,
    y[:,0] / n,
    label="cell lists",
)

plt.scatter(
    n,
    y[:,1] / n,
    label="cell lists, N3L",
)

n = 138 * x
plt.scatter(
    n,
    y[:,2] / n,
    label="neighbor lists",
)

plt.scatter(
    n,
    y[:,3] / n,
    label="neighbor lists, N3L",
)

plt.title("Performance of various algorithms on CPU")
plt.xlabel("Particles")
plt.ylabel("us per particle")
plt.legend()
plt.yscale("log")
plt.xscale("log")

plt.savefig("cpu")
