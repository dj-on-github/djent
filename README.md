# djent
djent is a reimplementation of the Fourmilab/John Walker random number test program ent.

The improvements are:

* Multiple input file names can be provided at once. This works nicely with the CSV format output.
* -h works are well as -u to get the help information.
* A more aggressive compression test. The existing test will report 0% compression on serially correlated data
* Bit reversal detection for serially correlated data.


 
