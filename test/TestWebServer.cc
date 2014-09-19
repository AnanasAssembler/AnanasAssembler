#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"

#include "base/WebServer.h"
#include <unistd.h>
#include <stdio.h>

#include <stdlib.h>

#include <unistd.h>

/*
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sys/wait.h>
#include <signal.h>
#include <netdb.h>
*/


class ReqHandler : public IRequestHandler
{
public:
  ReqHandler() {
    m_counter = 0;
  }

  virtual void OnInit() {
    m_pWebServer->AddRouteThrough(".jpeg");
    m_pWebServer->AddRouteThrough(".png");
    m_pWebServer->AddRouteThrough(".jpg");
  }
  
  virtual int OnRequest(struct mg_connection *conn, const string & uri, const string & id) {
    cout << "OnRequest. URI: " << uri << " ID: " << id << endl;
    if (uri == "/handle_post_request") {
      m_counter++;
      cout << "Counting: " << m_counter << endl;
      if (m_counter < 5) {
	m_pWebServer->SendHTMLFile(conn, "smorgas_server_data/please_wait.html");
      } else {
	m_pWebServer->SendHTMLFile(conn, "smorgas_server_data/Results.html");
	m_counter = 0;
      }
      return MG_TRUE;
    } else {
      m_pWebServer->SendHTMLFile(conn, "smorgas_server_data/Tangerine_index.html", "id=blah");
      //m_pWebServer->SendHTMLFile(conn, "Tangerine_index.html", "");
      return MG_TRUE;
    }
  }
 
  virtual void OnAuth(struct mg_connection *conn) {

  }

  virtual void OnOther(struct mg_connection *conn) {

  }

private:
  int m_counter;
};


int main( int argc, char** argv )
{
  /*
  commandArg<string> fileCmmd("-i","input file");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the MGWebServer.");
  P.registerArg(fileCmmd);
  
  P.parse();
  
  string fileName = P.GetStringValueFor(fileCmmd);
  */
  ReqHandler handler;

  WebServer server(&handler);

  
   
  return 0;
}
