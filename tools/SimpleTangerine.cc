
//#include <string>
//#include <stdio.h>
//#include <unistd.h>

#include "base/CommandLineParser.h"
//#include "util/SComm.h"
#include "base/FileParser.h"

int main(int argc,char** argv)
{


  for (int i=0; i<100000; i++) {
      
    char tmp[8000];
      
    cout << "Please paste in your sequence: " << endl;
    cin >> tmp;
    FILE * p = fopen("query.tangerine", "w");
    fprintf(p, ">YourSequence\n");
    fprintf(p, "%s\n", tmp);    
    
    fclose(p);
  }
  return 0;
}
