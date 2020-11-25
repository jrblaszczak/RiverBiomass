# River Productivity Modeling

## 2020-11-16 update

1. Issues with PM1 and PM2 were resolved with help from Bob S. (but I do want to claim the idea of using a rayleigh distribution for c which I think has helped constrain it considerably).

2. Where I need help from Bob H. is taking a second look at PM3 and the priors.



## 2020-11-13 Update

1. Currently working from "./code/Productivity_Model_Simulations - Oregon Test.Rmd" to simulate data

	* This RMD links to to the data simulation code and Stan models described below in (2)

2. There are three versions of the model to simulate data that correspond to Stan versions of those models in the code folder:

	* "Simulated_ProductivityModel1_Autoregressive.R" & "Stan_ProductivityModel1_Autoregressive.stan"

	* "Simulated_ProductivityModel2_Logistic.R" & "Stan_ProductivityModel2_Logistic.stan"

	* "Simulated_ProductivityModel3_ThinFilm.R" & "Stan_ProductivityModel3_ThinFilm.stan"

3. There are three versions of the data simulation R scripts which are in "*/code/Code to play around with params" which are there to play around with the input parameters to get a better handle on what's going on.

4. What I need help on from Bob:

	A. Please look at the "Simulated_ProductivityModel*" R scripts, specifically in relation to the process and observation error; you might create another folder for yourself in the code folder or open up the "Code to play around with params" folder to adjust the simulations

	B. Please take a look at the corresponding "Stan_ProductivityModel*" scripts and see if you notice errors, again in particular in relation to the process and observation error

	C. **In particular**, when looking at the data simulation for PM1, everytime I try to change the distribution for the process model to log normal (rlnorm) it goes insane (see "PM1_adjusting_params"), but if I leave it a normal distribution it's fine, but GPP is on the log scale so I'm pretty sure the process model should be log normally distributed but maybe I'm confused.


