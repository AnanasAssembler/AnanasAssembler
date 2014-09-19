#define FORCE_DEBUG

#include <string>
#include <stdio.h>
#include <unistd.h>

#include "base/CommandLineParser.h"
#include "util/SyncConn.h"
#include "base/FileParser.h"
#include "util/FindProcess.h"
#include "base/StringUtil.h"
#include "base/ThreadHandler.h"
#include <stdio.h>



class Notifier
{
public:
  Notifier() {
    m_counter = 0;
    m_pOut = NULL;
  }
  void ResetCounter() {
    m_counter = 0;
  }

  // Needs to handle sessions!!
  void SendResult(const string & result, const string & db) {
    m_mutex.Lock();
    if (result == "Server waiting") {
      //cout << "Caught server waiting..." << endl;
      //return;
    }
    cout << result << endl;
    if (m_pOut != NULL)
      fprintf(m_pOut, "%s", result.c_str());
    m_counter++;
    m_mutex.Unlock();
  }

  int Counter() {
    //m_mutex.Lock();
    return m_counter;
    //m_mutex.Unlock();
  }

  void OpenOut(const string & s) {
    //m_mutex.Lock();
    m_pOut = fopen(s.c_str(), "w");
    //m_mutex.Unlock();
  }
  void CloseOut() {
    //m_mutex.Lock();
    fclose(m_pOut);
    m_pOut = NULL;
    //m_mutex.Unlock();
  }

private:
  ThreadMutex m_mutex;
  int m_counter;
  FILE * m_pOut;

};


class DispatchThread : public IOneThread
{
public:
  DispatchThread(int port, const string & server, const string & name, Notifier * pNotifier) : m_ss(server, port)
  {
    m_port = port;
    m_name = name;
    m_pNotifier = pNotifier;
    cout << "Cond." << endl;
  }

  void Add(const string & s) {
    cout << "Adding." << endl;
    m_mutex.Lock();
    m_buffer.push_back(s);
    m_mutex.Unlock();
  }
 
protected:
  bool Pop(string & s) {
    s = "";
    if (m_buffer.isize() > 0) {
      s = m_buffer[0];
      for (int i=0; i<m_buffer.isize()-1; i++)
	m_buffer[i] = m_buffer[i+1];
      m_buffer.resize(m_buffer.isize()-1);
      return true;
    }
    return false;
  }

  virtual bool OnDie() {
     return true;
  }

  virtual bool OnDo(const string & msg) {
    cout << "Called OnDo " << msg << endl;
    Add(msg);
    /*string msg;
    while (true) {
      if (Pop(msg)) {
	string result;
	cout << "Sending." << endl;
	m_ss.SendRequest(result, tmp);
	cout << "Returned: " << endl;
	m_pNotifier->SendResult(result, m_name);
	//cout << result << endl;
      } else {
	usleep(100000);
      }
      }*/

    return true;
  }

  virtual bool OnInitialize(const string & msgExt) {
    cout << "Initializing." << endl;
    string msg;
    while (true) {
      if (Pop(msg)) {
	string result;
	cout << "Sending." << endl;
	m_ss.SendRequest(result, msg);
	cout << "Returned: " << endl;
	m_pNotifier->SendResult(result, m_name);
	//cout << result << endl;
      } else {
	usleep(100000);
      }
    }
    return true;
  }

private:
  string m_name;
  int m_port;
  ThreadMutex m_mutex;
  svec<string> m_buffer;
  Notifier * m_pNotifier;
  SyncConnClient m_ss;


};





class Dispatcher
{
public:
  Dispatcher() {
    m_inc = 4;
  }

  void ReadConfig(const string & s);
  void Set(const string & name, const string & db) {
    m_dbname.push_back(name);
    m_db.push_back(db);
  }

  void Feed(const string & msg) {
    cout << "Feed..." << endl;
    for (int i=0; i<m_threads.isize(); i++) {
      m_threads[i]->Add(msg);
    }
    /*
    for (int i=0; i<m_th.GetNumThreads(); i++) {
      m_th.Feed(i, msg);
      cout << "Fed." << endl;
      }*/
  }

  void StartServers(const string & exe, int baseport, const string & server) {
    int port = baseport;
    if (getProcessCount("ProtServer") == m_db.isize()) {
      cout << "Not starting servers, they are running already." << endl;
      m_port.resize(m_db.isize(), -1);
      for (int i=0; i<m_db.isize(); i++) {
	m_port[i] = port;
	port += m_inc;
      }
    } else {
      //int port = MYPORT;
      m_port.resize(m_db.isize(), -1);
      for (int i=0; i<m_db.isize(); i++) {
	string cmmd = exe + " -t " + m_db[i];
	cmmd += " -port " + Stringify(port) + " -server " + server + " > protserv.out." + Stringify(i) + " &";
	int r = system(cmmd.c_str());
	cout << "Starting server, please wait..." << endl;
	cout << cmmd << endl;
	m_port[i] = port;
	port += m_inc;
      }
    }
  }

  void StartThreads(Notifier * pNotifier, const string & server) {
    int i;
    // if (fork()) {
      for (i=0; i<m_db.isize(); i++) {
	cout << "Starting up thread." << endl;
	DispatchThread * p = new DispatchThread(m_port[i], server, m_dbname[i], pNotifier);
	m_threads.push_back(p);
	m_th.AddThread(p);
      }
      //for (i=0; i<m_db; i++) {
      //m_th.Feed(i, "do_it");
      //}
      //} else {
      //cout << "ERROR!!!" << endl;
      //}
      cout << "Threads started." << endl;
  }
  int NumThreads() const {return m_threads.isize();}
private:
  svec<string> m_dbname;
  svec<string> m_db;
  svec<int> m_port;
  int m_inc;
  ThreadHandler m_th;
  svec<DispatchThread*> m_threads; 

};

void Dispatcher::ReadConfig(const string & s)
{
  FlatFileParser parser;
  
  parser.Open(s);

  while (parser.ParseLine()) {
    if (parser.GetItemCount() < 2) 
      continue;
    m_dbname.push_back(parser.AsString(0));
    m_db.push_back(parser.AsString(1));
  }
 
}


int main(int argc,char** argv)
{

  commandArg<string> aStringCmmd("-server","name of the server running the host");
  commandArg<string> fileCmmd("-f","protein database fasta","/references/databases/transposon_lib/transposon_db.pep_noStop");
  commandArg<string> confCmmd("-c","config file", "");
  commandArg<string> exeCmmd("-e","server executable", "../ProtServer");
  commandArg<bool> interCmmd("-interactive","interactive, not server", false);

  commandLineParser P(argc,argv);
  P.SetDescription("Communication layer between protein alignment server(s) and the web server.");
  P.registerArg(aStringCmmd);
  P.registerArg(fileCmmd);
  P.registerArg(exeCmmd);
  P.registerArg(confCmmd);
  P.registerArg(interCmmd);

  P.parse();

  string server = P.GetStringValueFor(aStringCmmd);
  string fasta = P.GetStringValueFor(fileCmmd);
  string exe = P.GetStringValueFor(exeCmmd);
  string conf = P.GetStringValueFor(confCmmd);
  bool bInter = P.GetBoolValueFor(interCmmd);

  Dispatcher dd;
  if (conf != "")
    dd.ReadConfig(conf);
  else
    dd.Set("default", fasta);

  dd.StartServers(exe, MYPORT + 4, server);

  Notifier notifier;
   
  dd.StartThreads(&notifier, server);
  cout << "Threads are running." << endl;

    
  int i, j;

  if (bInter) {
    for (i=0; i<100000; i++) {
      
      char tmp[2048];
      
      cout << "Input sequence: " << endl;
      //string tmp;
      cin >> tmp;
      
      dd.Feed(tmp);
    }
    return 0;
  }

  SCommReceiver * pRec = GetReceiver(server.c_str());
  while (true) {
    bool bOK = false;
    char msg[4096];
    strcpy(msg, "");
    while (true) {
      if (pRec->Get(msg, sizeof(msg))) {
	bOK = true;
	cout << "Got message " << msg << endl;
	break;
      } else {
	//cout << "Waiting..." << endl;
	usleep(100000); 
      }
    }

    cout << "Loop exit." << endl;
    //if (!bOK)
    //continue;
     
    cout << "Parsing message" << endl;
    StringParser s;
    s.SetLine(msg);
    if (s.GetItemCount() != 2) {
      cout << "ERROR, Did not understand: " << msg << endl;
    }
    const string & input = s.AsString(0);
    const string & output = s.AsString(1);
    string done = output + ".done";
    
    stringstream tmp;
    FlatFileParser parser;
    //string tmptmp;
    cout << "Opening input " << input << endl;
    parser.Open(input);
    while (parser.ParseLine()) {
      tmp << parser.Line() << endl;
      //tmptmp += parser.Line();
      //tmptmp += "\n";
      //cout << parser.Line() << endl;
      //cout << tmptmp << endl;
    }

    //tmp << tmptmp;
    cout << "Done, sending." << endl;
    notifier.ResetCounter();
    notifier.OpenOut(output);
    cout << "Feeding" << endl;
    cout << tmp.str() << endl;
    dd.Feed(tmp.str());
    
    while (true) {
      if (notifier.Counter() > 0) {
	cout << "Counter=" << notifier.Counter() << endl;
      }
      if (notifier.Counter() == dd.NumThreads())
	break;
      usleep(100000);
    }
    cout << "Got it. " << endl;
    notifier.CloseOut();
    FILE * pDone = fopen(done.c_str(), "w");
    fprintf(pDone, "done\n");
    fclose(pDone);	
    

  }
  delete pRec;
 

  

  return 0;
}
