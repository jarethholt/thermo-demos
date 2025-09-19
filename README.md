# Thermodynamic systems demos

## Supercritical system

This demo comes from the Spectral Collective video [Analyzing a mean-field phase change](https://youtu.be/yEcysu5xZH0). It provides a tractable model of a supercritical phase transition.

For this system, energy is lower when the density is high, similar to the long-range attraction of molecules. If the volume is fixed and the system is in contact with a reservoir of both particles and energy, then the temperature $T$ and chemical potential $C$ are the primary thermodynamic variables. The free energy per unit volume is then

$$ F(d; T, C) = E - T S - C d \quad \text{where} \quad E(d) = -2 d^2, \quad S(d) = -d \log(d) - (1 - d) \log(1 - d) $$

as a function of the density $d$. The entropy term counts the number of distinct states available for the $d V$ particles.

The equilibrium density is the value $d_{eq}(T, C)$ that minimizes the free energy. Thus we need to find $d$ such that

$$ F_d = -4 d - C + T (\log(d) - \log(1 - d)) = 0, \quad F_{dd} = -4 + \frac{T}{d (1 - d)} > 0. $$

This problem has the natural range $d \in [0, 1]$. The equilibrium density is guaranteed to be $d_{eq} \in (0, 1)$ since

$$ F(d=0) = 0, \quad F_d(d=0) = \log(0^+) = -\infty $$
$$ F(d=1) = -2 - C, \quad F_d(d=1) = -\log(0^+) = +\infty. $$

That is, although the free energy is well-defined and finite at the boundaries, the values are approached with infinite derivatives that ensure any minima are in the interior.

Solving $F_{dd} > 0$ gives

$$ d < d_-(T) = \frac{1 - \sqrt{1 - T}}{2} < 1/2 \quad \text{or} \quad d > d_+(T) = \frac{1 + \sqrt{1 - T}}{2} > 1/2. $$

When $T > 1$ the function is concave on $(0, 1)$ and there is a unique minimum. This is the supercritical phase. For $T < 1$ the function is not concave. There will be a local minimum with $d < 1/2$, corresponding to a gas phase, if

$$ F_d(d_-) > 0 \quad \Rightarrow \quad C < C_{gas}(T) = -2 + 2 \sqrt{1 - T} + T \log\left( \frac{2 - T - 2 \sqrt{1 - T}}{T} \right). $$

Similarly, there will be a local minimum with $d > 1/2$ corresponding to a liquid phase if

$$ F_d(d_+) < 0 \quad \Rightarrow \quad C > C_{liq}(T) = -2 - 2 \sqrt{1 - T} + T \log\left( \frac{2 - T + 2 \sqrt{1 - T}}{T} \right). $$

There is a range of chemical potentials $C_{liq} < C < C_{gas}$ for which there are two local minima. One is the absolute minimum and the true equilibrium state; the other is a metastable state. The critical chemical potential separating liquid and gas equilibria turns out to be $C = -2$ at all temperatures.
