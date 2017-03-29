# djent
djent is a reimplementation of the Fourmilab/John Walker random number test program ent.

The improvements are:

* Multiple input file names can be provided at once. This works nicely with the CSV format output.
* -h works as well as -u to get the help information.
* The filename is present in CSV output
* The symbol size can be any number of bits up to 32. ent was constrained to 1 or 8.
* The SCC test can be either wrap-around or not wrap-around.
* The SCC result can be given a lag value to get a LAG-N correlation coefficient.
* A list of filenames to analyze can be read from a text file using -i filename.

Planned improvements are:

* Bit reversal options
* Adding a mode to give exactly the same output as ent for compatibilty with tools that parse the ent output.

```
./djent -h
Usage: djent [-b] [-l <n>] [-c] [-u] [-h] [-f] [-t] [filename] [filename2] ...

Compute statistics of random data.
  Author: David Johnston, dj@deadhat.com

  -l <n>    --symbol_length=<n> Treat incoming data symbols as bitlength n. Default is 8.
  -b        --binary            Treat incoming data as binary. Default bit length will be -l 1
  -c        --occurrence        Print symbol occurrence counts
  -w        --scc_wrap          Treat data as cyclical in SCC
  -n <n>    --lagn=<n>          Lag gap in SCC. Default=1
  -f        --fold              Fold uppercase letters to lower case
  -i <filename>  --inputfilelist=<filename> Read list of filenames from <filename>
  -t        --terse             Terse output
  -h or -u  --help              Print this text

```
  


 
