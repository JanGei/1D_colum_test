## 1-Dimensional Analytical Column Test

#### General Structure
This repository contains an interactive, analytical [model]([https://bokeh.org/](https://jangei.github.io/1D_colum_test_analytical/)) that solves 1-dimensional solute transport through a cylinder, representative of a column experiment. 
The basic structure of the model is developed in Python, while the interactive nature of this model functions through JavaScript. 
This is realized through the Python library [bokeh](https://bokeh.org/) which provides customizable JavaScript callbacks in the Python language.

#### Analytical Equation
![Analytical Equation](https://user-images.githubusercontent.com/99887101/195619550-84df8812-c116-453c-b32d-831872905382.PNG)

The analytical equation above is being solved through the model and contains the following elements:
- Normed concentration C(x,t)/C<sub>0</sub> of the solute 
- Location x within the column
- Time t
- First order rate constant λ
- Seepage velocity v
- Complimentary Error function erfc 
- Dispersion coefficient D
- Combined parameter H = 2λD / v<sup>2</sup>

#### User Interface

![AnalyticalUI](https://user-images.githubusercontent.com/99887101/195622686-1e3190a3-8ecf-486a-8a09-4605bf15db6a.PNG)

The UI of the model is essentially divided into four parts.  
The plot in the upper left part shows the normed solute concentration as a function of space, *i.e.* length of the column, at a given point in time.  
The plot in the lower left part shows a breakthroughcurve, *i.e.* the normed solute concentration as  a function of time at a given location. 
This location is specified by the grey diamond in the upper plot, which can be dragged in order to change its position.  
The upper right part of the UI contains the widgets through which the user can change the models mode, geometry, and the solutes properties.
Initially the user can choose between a continuous injection and a pulse injection. The picture above depicts the model in its 'Continuous Injection' state. 
If the mode of the model is switched to 'Pulse Injection', an additional slider appears, that controlls the duration of injection.
Furthermore, the user can choose between 'No Sorption' and 'Linear Sorption', as linear sorption of a solute can be modelled analytically. 
For this, the retardation factor R is computed:  
R  = 1 + (1-n)/n K<sub>d</sub> ρ<sub>d</sub>  
with:
- Porosity n
- Linear Partitioning Coefficient K<sub>d</sub>
- Solid Density ρ<sub>d</sub>

This retardataion factor is multiplied to the seepage velocity in order to model linear sorption. 
If the option 'Linear Sorption' is selected, three sliders, containing K<sub>d</sub> and ρ<sub>d</sub> appear.  
An extension of the problem to non-linear sorption can be found in this [repository]([https://bokeh.org/](https://github.com/JanGei/1D_column_test_numerical))

Below these options the 'Time' slider is found, through which the elapsed time within the simulation can be changed. 
Note that for practical reasons time is converted into pore volumes, but a conversion rate from pore volumes to hours is depicted in the slider title.
Following the 'Time' slider, all necessary and adjustable model parameters are displayed with a slider. 
Through the sliders the numerical values of the columns length and radius, the solutes first-order reaction coefficient, dispersion coefficient, and the setups flow rate and porosity can be altered.  

The lower right part of the UI contains a unit selection for the reaction, dispersion, and flow rate term. 
However, after selecting the unit of choice, the user needs to interact with the slider again, as this triggers the unit switch.  

Below the unit selection are two buttons that, when clicked, download a .csv file, containing x and y values of the upper and lower plot to the download folder of the computer.

#### Uncertainty

In experimental column setups the exact value of the first-order reaction coefficient and/or the dispersion coefficient is not known exactly.
To tackle this problem, the model allows the user to define a range with a lowest and highest possible value for both variables.
Within this range, a [latin hypercube](https://en.wikipedia.org/wiki/Latin_hypercube_sampling) is used to sample 100 parameter combinations containing both variables.
Thus, the model re-computes the aforeshown equation 100 times and computes the ensemble mean, as well as the upper and lower quartile and the minimum and maximum.
These values are also depicted in the figure in the upper left part of the UI.
Collapsing the range slider minimum and maximum into one point for both variables disables this option.

#### Future Work and Extensions

An extended (and numerical) version of this model can be found [here]([https://bokeh.org/](https://github.com/JanGei/1D_column_test_numerical)). 
It contains non-linear sorption (Freundlich and Langmuir). Further possible extensions and future work is discussed there.
