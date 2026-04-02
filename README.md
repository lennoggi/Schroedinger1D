# Schroedinger1D
A code evolving the 1D Schroedinger equation
$$
\begin{equation*}
    i\hbar\frac{\partial\Psi}{\partial t}\left(t, x\right) = H\left(t, x\right)\Psi\left(t, x\right)\equiv\left[-\frac{\hbar^2}{2m}\frac{\partial}{\partial x^2} + V\right]\Psi\left(t, x\right)
\end{equation*}
$$
with periodic boundary conditions and any real, time-independent potential $V\left(x\right)$ via the Lie-Trotter operator-splitting technique through a discrete Fourier transform. The code was built following [these lecture notes](https://virgilio.mib.infn.it/~oleari/public/laboratorio_fis_computaz/ris_num_eq_Schroedinger.pdf).
