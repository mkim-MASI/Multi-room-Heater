# Multi-room-Heater

Used MATLAB version 2022a

Executable is called "dynamics_generate"

No special instructions to run the execuatble, simply run it and the program will ask you for inputs (u = outside temp, hi, i={1, 2, 3, 4} are the initial heater placements). Make sure to have the sum of all hi be exactly 2 (ex: [0 1 1 0]). There are no checks for this, but you should do it for thr program to run properly. Also, choose a reasonable value for u (e.g. u = 6 is normal). Anything too large will result in the room temperatures becoming unbounded.

The file called "dynamics.slx" is our multi-room heating system with the scopes. Feel free to look at that in MATLAB to change values/parameters and see what happens.

"dynamics_generate.slx" is the file we used to generate our C code. The executable was created by first altering our generated C files to prompt the user for input and give output instead of using hard coded values. Then, we built the executable from these files.

'FinalReport.docx' is our final report for this project
