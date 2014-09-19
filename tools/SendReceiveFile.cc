
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
  commandArg<string>  fileCmmd("-f","file name");

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
  SCommTransmitter * pTrans = NULL;
 
  int port1 = 3491;
  int port2 = 3493;

  if (aBool) {
    int n = 0;
    //cout << "Will be the server (ctrl+c to exit, or will exit after 10 sends)!" << endl;
 
    if (fileName != "") {
      FlatFileParser parser;
      pTrans = GetTransmitter(port1);
  
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
      //cout << "Sending all data, bytes=" << k << endl;
      pTrans->SendWait(data);
      delete pTrans;
      pTrans = NULL;
      delete [] data;
      //return 0;
    } else {
    
      char tmp[2048];
      
      //sprintf(tmp, "This is message number %d", n);
      n++;
      
      //cout << "Type message." << endl;
      //string tmp;
      cin >> tmp;
      
      //cout << "Waiting for message to be retrieved..." << endl;
      pTrans = GetTransmitter(port1);
 
      pTrans->SendWait(tmp);
      //cout << "it's gone!" << endl;
      delete pTrans;
      pTrans = NULL;
    }
    //delete pTrans;
    //cout << "Limit reached, exiting" << endl;


    SCommReceiver * pResponse = GetReceiver(aString.c_str(), port2);
    char * tmp = new char[max];
    for (i=0; i<5000; i++) {
      cout << "Listening" << endl;
      if (pResponse->Get(tmp, max)) {
	cout << tmp << endl;
	break;
      } else {
	cout << "Waiting." << endl;
	sleep(1);
      }
    }
    delete [] tmp;
    delete pResponse;
    return 0;
  }

  if (aString == "") {
    cout << "You must specify the host name!" << endl;
    return -1;
  }


  char * tmp = new char[max];
  
  SCommReceiver * pRec = GetReceiver(aString.c_str(), port1);

  for (i=0; i<50000; i++) {
    if (pRec->Get(tmp, max)) {
      //delete pRec;
      cout << tmp << endl;
      //delete [] tmp;
 
       
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
      //cout << "Sleeping." << endl;
      //sleep(2);
      cout << "Sending data back, bytes=" << k << endl;
      pTrans = GetTransmitter(port2);
      pTrans->SendWait(data);
      cout << "Done!!" << endl;
      delete pTrans;
      pTrans = NULL;
      delete [] data;

      delete [] tmp;
      delete pRec;
   //pRec = GetReceiver(aString.c_str(), port1);
      return 0;

    } else {
      //cout << "Receive error!" << endl;
      sleep(1);
    }
  }

  delete [] tmp;
  delete pRec;
}
