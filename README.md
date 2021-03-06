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
* Test condition details (Volts, temp, id etc.) can be parsed from the filename and included in output. 
* MCV Min Entropy is estimated in addition to Shannon Entropy. The symbol and entropy are both reported
* The longest run and the symbol in the longest run are reported. For 1 bit-per-symbol analysis, a p-value is computed of the probability of a uniform random bit sequence having a longest run length equal to or less than the meaured run length. 

```
djent -h
Usage: djent [-brRpcCuhds] [-l <n>] [-i <input file list filename>] [filename] [filename2] ...

Compute statistics of random data.
  Author: David Johnston, dj@deadhat.com

  -i <filename>  --inputfilelist=<filename> Read list of filenames from <filename>
  -p             --parse_filename           Extract CID, Process, Voltage and Temperature from filename.
                                            The values will be included in the output.
  -l <n>         --symbol_length=<n>        Treat incoming data symbols as bitlength n. Default is 8.
  -b             --binary                   Treat incoming data as binary. Default bit length will be -l 1
  -r             --byte_reverse             Reverse the bit order in incoming bytes
  -R             --word_reverse             Reverse the byte order in incoming 4 byte words
  -c             --occurrence               Print symbol occurrence counts
  -C             --longest                  Print symbol longest run counts
  -w             --scc_wrap                 Treat data as cyclical in SCC
  -n <n>         --lagn=<n>                 Lag gap in SCC. Default=1
  -f             --fold                     Fold uppercase letters to lower case
  -t             --terse                    Terse output
  -e             --ent_exact                Exactly match output format of ent
  -s             --suppress_header          Suppress the header in terse output
  -h or -u       --help                     Print this text

 Notes
   * By default djent is in hex mode where it reads ascii hex data and converts it to binary to analyze.
     In hex mode, the symbol length defaults to 8, so normal hex files can be treated as a representation
     of bytes. The symbol length can be changed to any value between 1 and 32 bits using the -l <n> option.
   * With the -b option djent switches to binary reads in each byte as binary with a symbol length of 1.
   * To analyze ascii text instead of hex ascii, you need djent to treat each byte as a separate symbol, so
     use binary mode with a symbol length of 8. I.E. djent -b -l 8 <filename>
   * By default djent treats the MSB of each byte as the first. This can be switched so that djent treats
     the LSB as the first bit in each byte using the -r option.
   * Terse output is requested using -t. This outputs in CSV format. The first line is the header. If
     multiple files are provided, there will be one line of CSV output per file in addition to the header.
     The CSV header can be suppressed with -s.
   * To analyze multiple files, just give multiple file names on the command line. To read data in from
     the command line, don't provide a filename and pipe the data in. <datasource> | djent
   * The parse filename option =p picks takes four patterns from the filename to include in the output,
     This is so that it is easy to plot test conditions that are commonly encoded in a filename.
     Fields are delimited by uderscores. The four patters for CID, process, Voltage and Temperature are:
     _CID-<componentID>_ , _PROC-<process info>_, _<x>p<y>V_ and _<x>p<y>C_ . 'p' is the decimal point.
   * To compute the statistics, djent builds a frequency table of the symbols. This can be displayed
     using the -c option. The size of this table is what limits the the maximum symbol size. For each
     of the 2^n symbols, a 64 bit entry in a table is created. So for n=32, that's 32GBytes so the ability
     to handle large symbol sizes is limited by the available memory and the per process allocation limit.
   * The serial correlation coefficient is not wrap around by default, meaning that it does not compare
     the last value in the data with the first. To get wrap around behaviour, use the -w option.
   * The Lag-N correlation coefficient can be computed by using the -n <n> option. This causes the SCC
     computation to compare each Xth symbol with the (X+n)th symbol instead of the (X+1)th symbol.
     If you use wrap around with Lag-N, then the wrap around will reach n bits further into the start
     of the sequence.
   * The byte reverse option -r reverses the order of bits within each byte. The word reverse option -R
     reverses the order of bytes within each 32 bit word, from 3,2,1,0 to 0,1,2,3. Both -R and -r can
     be used together. Using -R with a data that isn't a multiple of 32 bits long will get padded with
     zeros, which may not be what you want. A padding warning will be sent to STDERR.
   * Instead of providing data file names on the command line, djent can be told to read a list of files
     from a text file. The file must have one filename per line. Lines beginning with # will be ignored.
     Use the -i <filename> option to request that djent reads the file list from <filename>.

 Examples
   Print this help
     djent -h

   Analyze hex file from stdin
     cat datafile.hex | djent

   Analyze binary file
     djent -b datafile.bin

   Analyze several files with CSV output
     djent -t data1.hex data2.hex data3.hex

   Analyze ascii symbols - Read in binary and set symbol size to 8.
     djent -b -l 8  textfile.txt

   Analyze binary file with parsable filename.
     djent -b -t -p  rawdata_CID-X23_PROC-TTFT_1p2V_25p0C_.bin
```
  


 
