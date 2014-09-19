
#include <string>
#include <stdio.h>
#include <unistd.h>

#include "base/CommandLineParser.h"
#include "util/SComm.h"
#include "base/FileParser.h"

int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-hostname","name of the server running the host", "");
  commandArg<bool>  aBoolCmmd("-s","be the server",false);
  commandArg<string>  fileCmmd("-f","file name","");

  commandLineParser P(argc,argv);
  P.SetDescription("Testing the socket communication between computers.");
  P.registerArg(aStringCmmd);
  P.registerArg(aBoolCmmd);
  P.registerArg(fileCmmd);

  P.parse();

  string aString = P.GetStringValueFor(aStringCmmd);
  string fileName = P.GetStringValueFor(fileCmmd);
  bool aBool = P.GetBoolValueFor(aBoolCmmd);

  int i;
  int max = 10000000;
 
  if (aBool) {
    int n = 0;
    cout << "Will be the server (ctrl+c to exit, or will exit after 10 sends)!" << endl;
    SCommTransmitter * pTrans = NULL;
 
    if (fileName != "") {
      FlatFileParser parser;
      pTrans = GetTransmitter();
  
      parser.Open(fileName);
      char * data = new char[max];
      parser.Open(fileName);
      int k = 0;
      while (parser.ParseLine()) {
	const char * p = parser.Line().c_str();
	int n = strlen(p);
	for (int x=0; x<n; x++) {
	  data[k] = p[x];
	  k++;
	}
	data[k] = '\n';
	k++;
      }
      data[k] = 0;
      cout << "Sending all data, bytes=" << k << endl;
      pTrans->SendWait(data);
      delete pTrans;
      pTrans = NULL;
      delete [] data;
      return 0;
    }
    


    for (i=0; i<100000; i++) {
      
      char tmp[256];
      
      sprintf(tmp, "This is message number %d", n);
      n++;
      
      cout << "Type message." << endl;
      //string tmp;
      cin >> tmp;
      
      cout << "Waiting for message to be retrieved..." << endl;
      pTrans = GetTransmitter();
 
      pTrans->SendWait(tmp);
      cout << "it's gone!" << endl;
      delete pTrans;
    }
    delete pTrans;
    cout << "Limit reached, exiting" << endl;
    return 0;
  }

  if (aString == "") {
    cout << "You must specify the host name!" << endl;
    return -1;
  }
  SCommReceiver * pRec = GetReceiver(aString.c_str());

  char * tmp = new char[max];

 
  for (i=0; i<50; i++) {
    if (pRec->Get(tmp, max))
      cout << "Got message: " << tmp << endl;
    else {
      cout << "Receive error!" << endl;
      sleep(1);
    }
  }

  delete [] tmp;
  delete pRec;
}
