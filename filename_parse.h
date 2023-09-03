/* look for vpattern in str. Return the match to found. Return True if found */

extern double voltage;
extern double temperature;
extern unsigned char deviceid[256];
extern unsigned char process[256];

int find_vpattern(char *str,char *found);
int find_tpattern(char *str,char *found) ;
int find_cidpattern(char *str,char *found); 
int find_procpattern(char *str,char *found); 
void parse_the_filename(char *filename);


