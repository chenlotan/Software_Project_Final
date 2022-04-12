# Software_Project_Final

changes made:
* the input of Jacobi Algorithm is now symmetric matrix (instead of running lnorm, gets lnorm matrix)
* the heuristic eigenvalues function was changed, some functions were added to sort the eigen. 
* the code for the module extension was written - build and import suppose to run (with the setup file)
* main function for the C code was added
* supporting csv and txt files in C is made in the same way - our code suppose to work for both of them
* etc

To do:
spkmeans.c:
1. spk not working - inaccurate numbers in Jacobi Algorithm and heuristic eigenvalues - !!! if k=1 in the end of this algorithm Error of type 2!!!
2. need to add validation that the matrix we get for Jacobi Algorithm is symmetric (?)
3. need to write comp.sh file

------------------------------
example for Jacobi matrix:
{{1.9997,1.2698,0.7305,0.0000},
{0.7071,-0.0065,-0.0116,0.7070},
{0.0009,0.7070,0.7070,0.0172},
{-0.0009,-0.7072,0.7070,0.0060},
{-0.7071,-0.0048,-0.0116,0.7070}}