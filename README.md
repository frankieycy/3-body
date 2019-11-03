## 3-Body Simulation

* Last: 2/11/2019
* 3-body simulation with simple Euler/predictor-corrector method
* a brief comparison of how the two methods err

## Code

* assign init conditions in `init()`
* easily extendable to a n-body sim
* change force field for a diff dynamical sim, e.g. molecular dynamics
* global params are used throughout, e.g. pos, vel, acc

## Outputs

* `pos.xyz`: coord of bodies
* `data.csv`: energies and errors over time
* `ke.png` & `pe.png`: comparison of energies under the two methods

## Discussion

* simple Euler: 1st-order scheme, errors rapidly accumulate. Approx 20% err over 4 earth periods (very bad)
* predictor-corrector: modify Euler prediction, signif improvememt of accuracy

## Other Schemes

* Euler-Cromer (simply swapping two lines but _far_ outperforms simple Euler)
* RK4
* many more...