# djent
djent is a reimplementation of the Fourmilab/John Walker random number test program ent.

The improvements are:

* Multiple input file names can be provided at once. This works nicely with the CSV format output.
* -h works as well as -u to get the help information.
* The filename is present in CSV output
* The symbol size can be any number of bits up to 32. ent was constrained to 1 or 8.
* The SCC test can be either wrap-around or not wrap-around.

Planned improvements are:

* Bit reversal options
* Adding Lag-N Correlation Coefficient
* Adding a mode to give exactly the same output as ent for compatibilty with tools that parse the ent output.

```
./djent -h
Usage: djent [-b] [-l <n>] [-c] [-u] [-h] [-f] [-t] [filename] [filename2] ...

Compute statistics of random data.
  Author: David Johnston, dj@deadhat.com
  -l <n>    --symbol_length=<n> Treat incoming data symbols as bitlength n. Default is 8.
  -b        --binary            Treat incoming data as binary. Sets default symbol length to 1
  -c        --occurence         Print symbol occurence counts
  -w        --scc_wrap          Treat data as cyclical in SCC
  -f        --fold              Fold uppercase letters to lower case
  -t        --terse             Terse output
  -h or -u  --help              Print this text
```
  


 
