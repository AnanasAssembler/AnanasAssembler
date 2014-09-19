#include "base/WebServer.h"
#include "base/FileParser.h"

WebServer * pWeb = NULL;

void Convert(string & o, int i) {
  char tmp[256];
  sprintf(tmp, "%d", i);
  o = tmp;
}

static int ev_handler(struct mg_connection *conn, enum mg_event ev) {
  return pWeb->EventHandler(conn, ev);
}

void IRequestHandler::OnRouteThrough(struct mg_connection *conn, const string & uri) {
  m_pWebServer->SendBinaryFile(conn, uri);
}

WebServer::WebServer(IRequestHandler * pHandler, int port)
{
  pWeb = this;
  m_pHandler = pHandler;
  m_pHandler->SetServer(this);
  m_pHandler->OnInit();
  server = mg_create_server(NULL, ev_handler);

  string p;
  Convert(p, port);
  mg_set_option(server, "listening_port", p.c_str());

  // Serve request. Hit Ctrl-C to terminate the program
  printf("Starting on port %s\n", mg_get_option(server, "listening_port"));
  for (;;) {
    mg_poll_server(server, 1000);
  }
}

WebServer::~WebServer()
{
  mg_destroy_server(&server);
}

void Split(string & base, string & ext, const string & uri)
{
  StringParser p;
  p.SetLine(uri, "?");
  base = p.AsString(0);
  ext = "";
  if (p.GetItemCount() > 1)
    ext = p.AsString(1);
  //cout << "Split: " << uri << " -> " << base << " " << ext << endl;
}

void Insert(string & line, const string & post_id)
{
  StringParser p;
  p.SetLine(line);
  int i;
  string n;

  if (p.GetItemCount() == 0)
    return;
  n = p.AsString(0);
  bool b = false;
  for (i=1; i<p.GetItemCount(); i++) {
    n += " ";
    if (p.AsString(i-1) == "method=\"POST\"" &&
	strstr(p.AsString(i).c_str(), "action=") != NULL) {
      char tmp1[256];
      strcpy(tmp1, p.AsString(i).c_str());
      tmp1[strlen(tmp1)-2] = 0;
      string tmp = tmp1;
      tmp += "?";
      tmp += post_id;
      tmp += "\">";
      n += tmp;
      b = true;
    } else {
      n += p.AsString(i);
    }
  }
  if (b) {
    //cout << "Replacing " << line << " with " << n << endl;
    line = n;
  }
}


int WebServer::EventHandler(struct mg_connection *conn, enum mg_event ev) {
  
  if (ev == MG_REQUEST) {
    string uri = conn->uri;
    int i;
    bool bRoute = false;
    for (i=0; i<m_through.isize(); i++) {
      if (strstr(uri.c_str(), m_through[i].c_str()) != NULL) {
	bRoute = true;
	break;
      }
    }
    if (bRoute) {
      if (strlen(uri.c_str()) > 0 && conn->uri[0] == '/')
	uri = &conn->uri[1];
      m_pHandler->OnRouteThrough(conn, uri);
    } else {
      //string base, ext;
      //Split(base, ext, uri);
      string ext;
      if (conn->query_string != NULL)
	ext = conn->query_string;
      return m_pHandler->OnRequest(conn, uri, ext);
    }
    return MG_TRUE;
  } else if (ev == MG_AUTH) {
    //printf("Auth\n");
    m_pHandler->OnAuth(conn);
    return MG_TRUE;
  } else {
    m_pHandler->OnOther(conn);
    //printf("Other\n");
    return MG_FALSE;
  }
}



void WebServer::AddRouteThrough(const string & s)
{
  m_through.push_back(s);
}

bool WebServer::SendData(struct mg_connection *conn, char * data, int len)
{
  mg_send_data(conn, data, len);
  return true;
}

bool WebServer::SendData(struct mg_connection *conn, const string & s)
{
  //cout << s << endl;
  mg_send_data(conn, s.c_str(), 1+strlen(s.c_str()));
  return true;
}

bool WebServer::SendBinaryFile(struct mg_connection *conn, const string & fileName)
{
  int data_len = 10000000;
  char * data = new char[data_len];
  
  FILE * pImage = fopen(fileName.c_str(), "rb");
  if (pImage == NULL) {
    cout << "ERROR: File not found: " << fileName << endl;
    return false;
  }
  int rr2 = fread(data, 1, data_len, pImage);	
  data_len = rr2;
  // printf("Read: %d\n", rr2);
  fclose(pImage);
  mg_send_data(conn, data, data_len);
  delete [] data;

  return true;
}



bool WebServer::SendHTMLFile(struct mg_connection *conn, const string & fileName, const string & post_id)
{
  FlatFileParser parser;
  
  parser.Open(fileName);
  int i;
  while (parser.ParseLine()) {
    string line = parser.Line();
    if (post_id != "")
      Insert(line, post_id);
    SendData(conn, line);
  }
  return true;
}

bool WebServer::GetVariable(struct mg_connection *conn, string & out, const string & name)
{
  // HARD CODED!!!!
  int max = 100000;
  char * buffer = new char[100000];
  buffer[max-1] = 0;

  int x = mg_get_var(conn, name.c_str(),
		     buffer, max-1);
  out = buffer;

  return true;
}

 

