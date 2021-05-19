------------------------------------------------------------------------------------------------------------------------------------------------------
--------------------------- REPRODUCING CAVITATION AREA USING PENALTY OPERATOR AND VARIATIONAL FORMULATION -------------------------------------------
------------------------------------------------------------------------------------------------------------------------------------------------------

FEM code written to reproduce cavitation area that generates in industrial application connected to bearings. Indeed in contact areas between bearings and the inner
pipe lubricant oil is subject to such high pressure that it undergoes cavitation. This small software is able to predict how big is the area in which cavitation occurs.


LABO.m -> main file

---------------------------------------------------------------------------PART 1-----------------------------------------------------------------------
Launching it, reproduce results of numerical test for selected variable:
			- Test1 reference case
			- Test2 problem with exact solution

for refinment values of 3,4,5 and 6

---->bool variable activation:
				- debug_plot shows penalty operator in action 
				- incr_debug check incremental error
				- cavitat_area_debug shows evolution of cavitation area



---------------------------------------------------------------------------PART 2-----------------------------------------------------------------------

Convergence test with epsilon = h^2


---------------------------------------------------------------------------PART 4-----------------------------------------------------------------------

Material EXTRA on Lagrangian Augmented (useful for quadratic BC or for stabilizing problems whose egenvectors were very close to each other)


--->Folder Numerical Error: Convergence test results

--------> SEE REPORT FOR MORE DETAILS


