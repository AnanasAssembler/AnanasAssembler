#ifndef WEBSERVER_H_
#define WEBSERVER_H_

#include "extern/mongoose/mongoose.h"
#include <string>
#include "base/SVector.h"

class WebServer;


inline string IDString(struct mg_connection *conn, const string & id) {
  string result = "ID";
  result += conn->remote_ip;
  result += "-";
  result += id;
  //result += conn->local_ip;
  return result;
}

class IDSystem
{
public:
  IDSystem() {}

  void Add(const string & id) {
    if (Index(id) != -1)
      return;
    m_ids.push_back(id);
    m_counter.push_back(0);
    cout << "Added ID counter: " << id << endl;
  }
  void Inc(const string & id) {
    int i = Index(id);
    if (i < 0) {
      cout << "ERROR: ID not found: " << id << endl;
      return;
    }
    m_counter[i]++;
  }

  int Counter(const string & id) {
    int i = Index(id);
    if (i < 0) {
      Add(id);
      //cout << "ERROR: ID not found: " << id << endl;
      return 0;
    }
    return m_counter[i];
  }

  void Remove(const string & id) {
    int i = Index(id);
    if (i < 0) {
      cout << "ERROR: ID not found: " << id << endl;
      return;
    }
    m_counter[i] = m_counter[m_counter.isize()-1];
    m_ids[i] = m_ids[m_ids.isize()-1];
    m_counter.resize(m_counter.isize()-1);
    m_ids.resize(m_ids.isize()-1);
    cout << "REMOVED ID counter: " << id << endl;
  }

private:
  int Index(const string & id) const {
    int i;
    for (i=0; i<m_ids.isize(); i++) {
      if (m_ids[i] == id)
	return i;
    }
    return -1;
  }
  svec<int> m_counter;
  svec<string> m_ids;
};




// Currently, this is a very thin layer...
class IRequestHandler
{
 public:
  IRequestHandler() {m_pWebServer = NULL;}
  virtual ~IRequestHandler() {}

  virtual void OnInit() = 0;
  virtual int OnRequest(struct mg_connection *conn, const string & uri, const string & id) = 0;
  virtual void OnAuth(struct mg_connection *conn) = 0;
  virtual void OnOther(struct mg_connection *conn) = 0;

  virtual void OnRouteThrough(struct mg_connection *conn, const string & uri);
 
  void SetServer(WebServer * p) {
    m_pWebServer = p;
  }
 protected:
  WebServer * m_pWebServer;
};


class WebServer
{
 public:
  WebServer(IRequestHandler * pHandler, int port = 8080);
  ~WebServer();
  
  bool SendData(struct mg_connection *conn, char * data, int len);
  bool SendData(struct mg_connection *conn, const string & s);
  bool SendBinaryFile(struct mg_connection *conn, const string & fileName);
  bool SendHTMLFile(struct mg_connection *conn, const string & fileName, const string & post_id = "");

  bool GetVariable(struct mg_connection *conn, string & out, const string & name);

  // Routes through objects (useful for putting up images etc.)
  void AddRouteThrough(const string & s);

  int EventHandler(struct mg_connection *conn, enum mg_event ev);

 private:
  struct mg_server *server;
  IRequestHandler * m_pHandler;
  svec<string> m_through;
};











#endif //WEBSERVER_H_

