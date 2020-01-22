# Abaqus-VUMAT-Johnson_Cook

This VUMAT subroutine is to implement the Johnson-Cook model for explicit systems. JC model encapsulates the dependency of plastic strain, plastic strain rate and temperature.


The program contains iteration based on the the reduction of boundary instead of commonly used Newton-Raphson method. This modification makes the program robust at the expense of running time. The plastic strain history is used to make the program efficient.
