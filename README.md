# Phage Listeria Analysis Model
This code is used to model the effects of temperature and pH on *Listeria monocytogenes'* sensitivity to bacteriophage (phage) treatment on a lab-scale Hispanic-style fresh cheese. 

## How to use

Raw data: Cheese_Phage_Masterfile.csv
Code: Phage_Listeria_Analysis_Model.R

We constructed two linear mixed effects models for **temperature** and **pH** using the "lmer" function in the "lme4" R package. For each model:

Response: log of the number (log_count) of *L. monocytogenes*

Random effects: 
- replicates (rep)
- plate nested within milk batch

Fixed effects:
- temperature or pH
- day
- phage
- strain
- age of the milk (milk_age)
- log of the aerobic plate counts (bacterial counts in the milk before cheese was made; milk_apc)

## How to Cite
Henderson, L.O.,Cabrera-Villamizar, L.A., Skeens, J., Kent, D., Murphy, S,. Wiedmann, M., and Guariglia-Oropeza, V. 2019. Dairy-specific environmental conditions and serotype affect *Listeria monocytogenes* susceptibility to phage treatment in a Hispanic-style fresh cheese model. Journal of Dairy Science. *In submission*

## Authors
L. O. Henderson, L. A. Cabrera-Villamizar, J. Skeens, D. Kent, S. Murphy, M. Wiedmann, and V. Guariglia-Oropeza
