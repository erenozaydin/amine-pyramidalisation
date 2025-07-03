# amine-pyramidalisation
This program takes in a .xyz file (which contains atomic coordinates) for a molecule containing a singular nitrogen atom and calculates a value for PΣ (a measure of pyramidality about the N atom; the extent of deviation from planarity).
### INPUT
A .xyz file formatted like so: 
<pre>
  (number of atoms in the molecule) 
  (comment line) 
  (element symbol) 
  (x co-ordinate) 
  (y co-ordinate) 
  (z co-ordinate) 
  eg. N -2.58458874 2.08215282 0.01332863 </pre>

The file must contain **exactly one nitrogen atom** bonded to at least three other atoms. 
### OUTPUT
The program prints the calculated PΣ value to 10 significant figures to the terminal.
### EXAMPLE OUTPUT

<pre> Please enter the full path for the .xyz file: C:\Users\YourName\filename.xyz

P_Sigma = (number to 10 significant figures)

Would you like to go again with another .xyz file? (yes/no): no
Thank you, I hope you enjoyed this work - Eren Ozaydin :) <pre>
