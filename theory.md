# Finite Difference Methods for Options Pricing

The Black-Scholes PDE is

$$
\boxed{
\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^{2}S^{2}\frac{\partial^{2}V}{\partial S^{2}}+(r-D)S\frac{\partial V}{\partial S} - rV = 0}
$$

where $V=V(S,t)$ is the price of an option on a stock that pays a continuous dividend yield $D$, has a strike $K$, expires $T$ years from now, and changes with stock price $S$ and time $t$. This equation has an analytical solution, however we can also solve for V numerically by approximating the derivative terms using finite differences. For example, we can approximate the time derivative as

$$
\frac{\partial V}{\partial t} = \lim_{\delta t \to 0}\frac{V(S,t+\delta t) - V(S,t)}{\delta t}
$$

To implement the finite differences computationally, we first choose $S^{*} >> K$ and define a grid by dividing the range $[0,S^{*}]$ into $N$ price steps and the range $[0,T]$ into $M$ time steps such that their lengths are

$$
\delta S = \frac{S^{*}}{N}
$$

and

$$
\delta t = \frac{T}{M}
$$

$S$ values and $t$ values such that

$$
\boxed{
V_{n}^{m}=V(n\delta S,m \delta t), \quad n \in [0,1,2...,N], \quad m \in [0,1,2,...,M]}
$$

So, if we are backward marching, we have

$$
\begin{align*}
\frac{\partial V}{\partial t} &= \lim_{\delta t \to 0}\frac{V(S,t+\delta t) - V(S,t)}{\delta t}
\end{align*}
$$

$$
\boxed{
    \begin{aligned}
\frac{\partial V}{\partial t}\approx\frac{V_{n}^{m}-V_{n}^{m-1}}{\delta t}
\end{aligned}
}
$$

Note here that we could define $V_{n}^{m}$ as

$$
V_{n}^{m}=V(n\delta S,T-m \delta t)
$$

in which case we would be forward marching. In this example, we will be backward marching.

Now lets approximate the derivatives with respect to $S$ using central finite differences. To do this we will use the taylor expansions of $V$, which are

$$
V(S+\delta S, t) = V(S,t)+ \frac{\partial V}{\partial S}\delta S+\frac{1}{2}\frac{\partial^{2} V}{\partial S^{2}} \delta S^{2} + O(\delta S^{3})
$$

and

$$
V(S-\delta S, t) =V(S,t)- \frac{\partial V}{\partial S}\delta S+\frac{1}{2}\frac{\partial^{2} V}{\partial S^{2}} \delta S^{2} + O(\delta S^{3})
$$

Subtracting these two expressions we get

$$
V(S+\delta S, t) - V(S-\delta S, t) \approx 2\frac{\partial V}{\partial S}\delta S\\
\frac{\partial V}{\partial S} \approx \frac{V(S+\delta S, t) - V(S-\delta S, t)}{2\delta S}
$$

Substituting our grid notation we get

$$
\boxed{
\frac{\partial V}{\partial S} \approx \frac{V_{n+1}^{m} - V_{n-1}^{m}}{2\delta S}}
$$

Now, if we add the two expressions together, we get

$$
V(S+\delta S, t) + V(S-\delta S, t) \approx 2V(S,t)+\frac{\partial^{2} V}{\partial S^{2}}\delta S^{2}
$$

$$
\frac{\partial^{2} V}{\partial S^{2}}\approx\frac{V(S+\delta S, t)- 2V(S,t) + V(S-\delta S, t)}{\delta S^{2}}
$$

Substituting our grid notation we get

$$
\boxed{
\frac{\partial^{2} V}{\partial S^{2}}\approx\frac{V_{n+1}^{m}- 2V_{n}^{m} + V_{n-1}^{m}}{\delta S^{2}}}
$$

Now, we note that $S^{2} = n^{2}\delta S^{2}$ and $S=n \delta S$, and we substitute the approximations for the derivatives back into the Black-Scholes PDE to get

$$
\frac{V_{n}^{m}-V_{n}^{m-1}}{\delta t} + \frac{1}{2}\sigma^2 n^2(V_{n+1}^{m}- 2V_{n}^{m} + V_{n-1}^{m})+(r-D)n\frac{V_{n+1}^{m} - V_{n-1}^{m}}{2} -rV_{n}^{m}= 0
$$

Solving for $V_{n}^{m-1}/ \delta t$ we get

$$
\begin{align*}
\frac{V_{n}^{m-1}}{\delta t} = V_{n-1}^{m}\left[ \frac{1}{2}\left(\sigma^{2}n^{2}-(r-D)n\right)\right] +&\\
V_{n}^{m}\left[\frac{1}{\delta t}-\sigma^{2}n^{2}-r\right] +& \\V_{n+1}^{m}\left[\frac{1}{2}\left(\sigma^{2}n^{2}+(r-D)n\right)\right]
\end{align*}
$$

and so

$$
\boxed{
V_{n}^{m-1} = A_{n}V_{n-1}^{m} + B_{n}V_{n}^{m}+C_{n}V_{n+1}^{m}
}
$$

where

$$
\boxed{
\begin{align*}
A_{n} &=\frac{1}{2}\left(\sigma^{2}n^{2}-n(r-D)\right)\delta t\\
B_{n} &=1-\left(\sigma^{2}n^{2}+r\right)\delta t \\
C_{n} &= \frac{1}{2}\left(\sigma^{2}n^{2}+n(r-D)\right)\delta t
\end{align*}
}
$$

**Boundary Conditions**

At the boundaries of the grid, the Black-Scholes boundary conditions should hold, i.e. the following should be true:

1. At expiry, the value should be equal to the payoff, i.e.

$$
C(S,T) = V_{n}^{M}=\text{max}(S-K,0) = \text{max}(n \delta S-K,0)
$$

for call options, and

$$
P(S,T) = V_{n}^{M}=\text{max}(K-S,0) = \text{max}(K-n \delta S,0)
$$

for put options.

2. Given that the Black Scholes framework models the stock price with Geometric Brownian Motion, i.e.

$$
dS_{t} = \mu Sdt+\sigma Sd W_{t}
$$

whenever $S$ falls to 0 between $t$ and $T$, it will remain at 0. and so

$$
V(0,t) = 0
$$

3. Recall that the closed form value of a call option under the Black Scholes framework is given by

$$
V(S,t) = Se^{-D(T-t)}N(d_{1}) -Ke^{-r(T-t)}N(d_{2})
$$

where

$$
\begin{align*}
d_{1} &= \frac{\ln(S/K)+(r-D+\sigma^{2}/2)(T-t)}{\sigma \sqrt{T-t}}\\
d_{2}&=\frac{\ln(S/K)+(r-D-\sigma^{2}/2)(T-t)}{\sigma \sqrt{T-t}}\
\end{align*}
$$

For large $S$ we have

$$
N(d_{1})_{S \to \infty} = N(d_{2})_{S \to \infty} = 1
$$

and so we have

$$
V(S^{*},t) = S^{*}e^{-D(T-t)} -Ke^{-r(T-t)}
$$
