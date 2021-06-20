# Finite Element analysis of the Reynolds lubrication equation with cavitation


FEM code able to reproduce cavitation area generated in industrial applications in lubricatn oil between bearings and a rotor shaft: in contact areas between bearings
and the shaft, lubricant oil is subject to such high pressure that it undergoes cavitation. 
Using a penalty operator algorithm applied to the weak formulation of the simplified Navier-Stokes set of equations, the software is able to predict the size of the contact
area in which cavitation happens. Further Augmented Lagrangian Algorithm has been implemented to account for instabilities arising by the use of non-linear boundary conditions.

**Installation**\
\
Git pull the repository and be sure you have MATLAB2018 installed

**Execution**\
\

(PART 1)
Open and launch LABO.m  (main file)
Launching it, reproduce results of numerical test for selected variable: "Test1" is the reference case whereas "Test2" is problem with exact solution

bool variable activation:
				- debug_plot shows penalty operator in action 
				- incr_debug check incremental error
				- cavitat_area_debug shows evolution of cavitation area

(PART 2)

Convergence test with epsilon = h^2

(PART 4)

Material EXTRA on Lagrangian Augmented (useful for quadratic BC or for stabilizing problems whose egenvectors were very close to each other)


--->Folder Numerical Error: Convergence test results

--------> SEE REPORT FOR MORE DETAILS


