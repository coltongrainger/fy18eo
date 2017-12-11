### Prob 36

We'll show that the 2nd order non-homogenous DE
$$
x'' + 4x = f(t)
$$
has solution $x(t) = \frac12 \int_0^t f(t-\tau) \sin 2\tau \, d\tau$.

Applying the $\mathcal{L}$ operator to both sides of the DE, then isolating $X(s)$, we have

$$
X = \frac{F(s)}{s^2 + 4} = \frac12 \cdot F(s) \cdot \frac{2}{s^2 + 4}.
$$

By the convolution theorem,

$$
\frac12\mathcal{L}^{-1}\left\{F(s) \cdot 2/(s^2 + 4))\right\} 
= \frac12 \int_0^t f(t-\tau) \sin 2\tau \, d\tau$.
$$

QED.
