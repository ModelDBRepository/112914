# MATLAB code for multistable anti-phase bursting

The MATLAB scripts posted below reproduce the simulation results corresponding to Figs. 1, 10 and 11 of the manuscript:

| Victor Matveev, Amitabha Bose, Farzan Nadim
**Capturing the bursting dynamics of a two-cell inhibitory network using a one-dimensional map** (2007)
*Journal of Computational Neuroscience*, **23**: 169.
[Full Text](http://www.springerlink.com/content/t4v879k50364471u/) |

Please place all files in the same directory before running the main simulation script **MultiMovie.m**

---

|  |  |
|---|---|
| ♦ [MultiMovie.m](http://web.njit.edu/%7Ematveev/Burst/MultiMovie.m) | This m-script illustrates the geometry of bursting corresponding to Figs. 1, 10 and 11, in the form of a movie, showing both the time course and the phase-plane dynamics of the model. The program prompts the user to select the number of spikes per burst. Only the initial condition is affected by user input; all parameters are set to fixed values listed in Appendix A1.  |
| ♦ [burstODE.m](http://web.njit.edu/%7Ematveev/Burst/burstODE.m) | This ODE m-file implements the model equations (Eqs. 6). The file follows the standard format used by the built-in MATLAB ode integrators. |
| ♦ [Vnullcline.m](http://web.njit.edu/%7Ematveev/Burst/Vnullcline.m) | Calculates the V-nullcline; used by the main MultiMovie script above. |

---

Supported in part by the **National Science Foundation** grants
**DMS 0417416** (Victor Matveev), **DMS 0615168** (Amitabha Bose)
and the **National Institutes of Health** grant **MH-60605** (Farzan Nadim).

---

Victor Matveev
[http://web.njit.edu/~matveev](http://web.njit.edu/%7Ematveev)

This server is running a
[Redhat](http://www.redhat.com/) distribution of
[Linux](http://www.linux.org/).

Last modified: Feb 4, 2007

---

2025-06-02: Standardized to Markdown.