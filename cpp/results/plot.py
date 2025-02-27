import csvp
import matplotlib.pyplot as plt


hdr, x, y = csvp.parse("benchmark.out"," ",float)

y /= 18
y *= 1e6

x = (x ** 3)[:,0]

n = 80 * x

plt.scatter(
    n,
    y[:,0] / n,
    label="cell lists",
)

n = 138 * x
plt.scatter(
    n,
    y[:,1] / n,
    label="neighbor lists",
)

plt.title("Performance of various algorithms on CPU")
plt.xlabel("Particles")
plt.ylabel("uSeconds per Particle per Timestep")
plt.legend()
plt.xscale("log")

plt.savefig("cpu")
