#define FORCE_DEBUG

#include <string>
#include "base/CommandLineParser.h"
#include "base/FileParser.h"
#include "base/ThreadHandler.h"
#include <math.h>

class SharedData
{
public:
  SharedData() {
    m_chunk = 4000000;
    //m_data.resize(m_chunk);
    m_k = 0;
    //m_data.resize(1);
  }

  void Add(int val, int howmany = 1) {
    m_mutex.Lock();
    for (int i=0; i<howmany; i++)
      m_data.push_back(val);
    //if (m_k >= m_data.isize())
    //m_data.resize(m_data.isize() + m_chunk);
    //m_data[m_k] = val;
    //m_k++;
    m_mutex.Unlock();
  }

  int Size() const {return m_data.isize();}
  int Get(int i) const {
    return m_data[i];
  }

  void Print() {
    for (int i=0; i<m_data.isize(); i++) 
      cout << i << "\t" << m_data[i] << endl;
  }

private:
  svec<int> m_data;
  ThreadMutex m_mutex;
  int m_k;
  int m_chunk;
};


class MyThread : public IOneThread
{
public:
  MyThread() {
    m_pData = NULL;
    m_id = -1;
    m_n = 0;
    m_howmany = 1;
  }

  MyThread(SharedData * p, int id, int n) {
    m_pData = p;
    m_id = id;
    m_n = n;
    m_howmany = 10;
  }

  void SetData(SharedData * p, int id) {
    m_pData = p;
    m_id = id;
  }
protected:

  virtual bool OnDie() {
    cout << "Killed!!" << endl;
    return true;
  }

  virtual bool OnDo(const string & msg) {
    cout << "Doing stupid things w/ " << msg << endl;
    for (int i=0; i<m_n; i++) {
      if ((i + 1) % m_howmany == 0)
	m_pData->Add(m_id, m_howmany);
      //cout << m_pData->Size() << endl;
      double v = 0.;
      for (int j=0; j<250; j++) {
	double f = exp(-(double)i);
	f = exp(-f);
	f = exp(-f);
	v += f;
      }
    }
    cout << "Done!" << endl;
    return true;
  }

  virtual bool OnInitialize(const string & msg) {
    cout << "Initializing!" << endl;
    cout << "Done init!" << endl;
    return true;
  }
private:
  SharedData * m_pData;
  int m_id;
  int m_howmany;
  int m_n;
};



int main( int argc, char** argv )
{

  commandArg<int> nCmmd("-n","number of threads");
  commandLineParser P(argc,argv);
  P.SetDescription("Testing the thread handler.");
  P.registerArg(nCmmd);
  
  P.parse();
  
  int n = P.GetIntValueFor(nCmmd);
  
  /*MyThread one, two;
  one.Initialize("something");
  one.Do("ah well");
  one.Do("ah well, whatever");

  two.Initialize("anything");
  two.Do("what the...");
  two.Do("even more");
   
  while (!one.Done()) {
    usleep(10000);
  }
  one.Die();

  cout << "All done, testing ThreadHandler." << endl;
  */
  int i;
  ThreadHandler th;

  SharedData data;
  cout << &data << endl;
  //data.Print();

  int size = 4000000;

  for (i=0; i<n; i++) {
    char tmp[256];
    sprintf(tmp, "%d", i);
    string init = "init_";
    init += tmp;
    cout << init << endl;
    th.AddThread(new MyThread(&data, i+1, size/n), init);    
  }
  for (i=0; i<n; i++) {
    char tmp[256];
    sprintf(tmp, "%d", i);
    string init = "feed_";
    init += tmp;
    cout << init << endl;
    th.Feed(i, init);
  }

  while (!th.AllDone()) {
    usleep(10000);
  }
  cout << data.Size() << endl;

  return 0;
}



