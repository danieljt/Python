The "Vil vite" Center in bergen has an interesting problem for young students.
Say you have the digits 0 - 9, and you want to place the digits in
the following equation.

(___ * __) - (___ * __) = ?

Example:
324 * 67 - 109*85 = 20836

At the "vil vite "center in bergen, the students are to use each of the
digits once in the equation, and the goal is place all 10 digits so
that the equation gives 0 as an answer. The total number of
possible combinations are 3628800, meaning the studets
have a lot of combinations to choose from.

In this program, we utilize a brute force algorithm where we try to
find all the possible solutions to the problem, and write them to a file.
By using simple and known algrithms, we can find all solutions of the
above equation that equal the number a. (where of course "a" in
the "Vil vite" center problem is 0).

EXAMPLE USAGE:
from numbertest import *   # Import program
run = numbertest()         # Set up instance
run.write("filename")      # Run the program and write to file
run.findvalue(5)           # Find all combinations that give 5 and print
