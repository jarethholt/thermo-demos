# Thermodynamic systems demos

## Supercritical system

This demo comes from the Spectral Collective video [Analyzing a mean-field phase change](https://youtu.be/yEcysu5xZH0). It provides a tractable model of a supercritical phase transition.

For this system, energy is lower when the density is high, similar to the long-range attraction of molecules. If the volume is fixed and the system is in contact with a reservoir of both particles and energy, then the temperature $T$ and chemical potential $C$ are the primary thermodynamic variables. The free energy per unit volume is then

$$ F(d; T, C) = E - T S - C d \quad \text{where} \quad E(d) = -2 d^2, \quad S(d) = -d \log(d) - (1 - d) \log(1 - d) $$

as a function of the density $d$. The entropy term counts the number of distinct states available for the $d V$ particles.

The equilibrium density is the value $d_{eq}(T, C)$ that minimizes the free energy. Thus we need to find $d$ such that

$$ F_d = -4 d - C + T (\log(d) - \log(1 - d)) = 0, \quad F_{dd} = -4 + T \left( \frac{1}{d} + \frac{1}{1 - d} \right) > 0. $$

This problem has the natural range $d \in [0, 1]$. The equilibrium density is guaranteed to be $d_{eq} \in (0, 1)$ since

$$ F(d=0) = 0, \quad F_d(d=0) = \log(0^+) = -\infty $$
$$ F(d=1) = -2 - C, \quad F_d(d=1) = -\log(0^+) = +\infty. $$

That is, although the free energy is well-defined and finite at the boundaries, the values are approached with infinite derivatives that ensure any minima are in the interior.

This system has a critical point at $(T_0, C_0) = (1, -2)$. There are three regions of the phase diagram. First, when $T < T_0$, there are two relative minima: a "gas" density $<1/2$ and a "liquid" density $>1/2$. For $C < C_0$ the gas density is the absolute minimum; for $C > C_0$ the liquid density is the absolute minimum. Along the equilibrium line $C = C_0$ the two minima have the same free energy and the phases coexist. For $T > T_0$ there is only one relative minimum (which is also the absolute minimum) and thus no phase change. While this density is $<1/2$ for $C < C_0$ and vice versa, it is not referred to as a gas or liquid density because it varies smoothly with $C$.

### Approximations

To have good starting points for the optimization, consider the limits of low and high densities. For low density:

$$ F_d \approx -C + T \log(d) = 0 \quad \Rightarrow \quad d_{eq} \approx e^{C/T} $$

For high density, $d = 1 - x$ with $x \ll 1$:

$$ F_d \approx -4 - T \log(x) = 0 \quad \Rightarrow \quad d_{eq} = 1 - e^{-4/T}.  $$
