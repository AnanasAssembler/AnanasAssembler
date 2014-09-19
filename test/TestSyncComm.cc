
#include <string>
#include <stdio.h>
#include <unistd.h>

#include "base/CommandLineParser.h"
#include "util/SyncConn.h"
#include "base/FileParser.h"

int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-hostname","name of the server running the host");
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
 
  char tmp[2048];

  if (aBool) {
    for (int i=0; i< 20; i++) {
      SyncConnServer server(aString);
      string request;
      server.WaitForRequest(request);
      cout << "GOT Request: " << request << endl;
      cout << "Type answer: ";
      cin >> tmp;
      server.SendResult(tmp);
    }
  } else {
    for (int i=0; i< 20; i++) {    
      cout << "Type request: ";
      cin >> tmp;
      
      SyncConnClient client(aString);
      string result;
      client.SendRequest(result, tmp);
      cout << "GOT Answer: " << result << endl;
    }
  }
  return 0;
}
