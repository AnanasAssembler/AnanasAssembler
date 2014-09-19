#ifndef _HMMDECODE_H_
#define _HMMDECODE_H_

#include "base/SVector.h"
#include <string>
#include <iostream>

#define HIGH_SCORE 99999999999.

//===============================================
class HMMStateScore
{
 public:
  HMMStateScore() {
    m_back = -1;
    m_score = HIGH_SCORE;
  }

  bool Merge(double score, int back) {
    //if (score < 0.)
    //score = -score;
    if (score < m_score) {
      m_score = score;
      m_back = back;
      return true;
    }
    return false;
  }


  void Reset() {
    m_back = -1;
    m_score = HIGH_SCORE;
  }
 
  void AddScore(double d) {
    m_score += d;
  }
    

  double GetScore() const {return m_score;}
  int GetBack() const {return m_back;}


 private:
  int m_back;
  double m_score;
};



//===============================================
class HMMHistoryNode
{
 public:
  HMMHistoryNode() {
    m_exit = -1;
    m_score = HIGH_SCORE;
    int m_back = -1;
    m_word = -1;
  }


  void Set(double score, int ex, int back, int word) {
    m_exit = ex;
    m_score = score;
    m_back = back;
    m_word = word;
  }


  int GetExit() const {return m_exit;}
  double GetScore() const {return m_score;}
  int GetBack() const {return m_back;}
  int GetWord() const {return m_word;}

 private:
  int m_exit;
  double m_score;
  int m_back;
  int m_word;
};



// The "trace-back" array
class HMMHistory
{
 public:
  HMMHistory() {
    m_lastFrame = 0;
    m_lastFrameStart = 0;
    m_index = 0;
  }


  void PossiblySize(int presize) {
    if (m_data.isize() > 0)
      return;
    m_data.resize(presize);
    m_lastFrame = 0;
    m_lastFrameStart = 0;
    m_index = 0;
  }

  int size() const {return m_data.isize();}
  const HMMHistoryNode & operator[] (int i) const {return m_data[i];}
  HMMHistoryNode & operator[] (int i) {return m_data[i];}

  int GetLastFrameStart() const {return m_lastFrameStart;}
  int GetLastFrameEnd() const {return m_index;}

  int Add(int ex, int back, double score, int word) {
    if (m_data.isize() <= m_index) {
      m_data.resize(m_index + 1000000);
    }

    m_data[m_index].Set(score, ex, back, word); 
    if (ex != m_lastFrame) {
      m_lastFrame = ex;
      m_lastFrameStart = m_index;
    }
    m_index++;
    return m_index-1;
  }


 private:
  svec<HMMHistoryNode> m_data;
  int m_lastFrame;
  int m_lastFrameStart;
  int m_index;
};



class HMMTracebackNode
{
 public:
  HMMTracebackNode() {
    m_exit = 0;
    m_entry = 0;
    m_score = 0;
    m_index = -1;
  }


  void Set(int entry, int ex, double score, const string & name, int index) {
    m_entry = entry;
    m_exit = ex;
    m_score = score;
    m_name = name;
    m_index = index;
  }


  int From() const {return m_entry;}
  int To() const {return m_exit;}
  double Score() const {return m_score;}
  const string & Name() const {return m_name;}
  int Index() const {return m_index;}

 private:
  int m_entry;
  int m_exit;
  double m_score;
  string m_name;
  int m_index;
};

class HMMStateFrameScoreSource
{
 public:
  virtual ~HMMStateFrameScoreSource() {}

  virtual double GetScore(const string & stateName) const = 0;

};

//===========================================================
class HMMStateFrameScores : public HMMStateFrameScoreSource
{
 public:
  HMMStateFrameScores() {
    m_index = 0;
  } 
  virtual ~HMMStateFrameScores() {}

  void Reset() {    
    m_index = 0;
  }

  void Set(const string & stateName, double score) {
    if (m_index >= m_names.isize()) {
      m_names.resize(m_index + 1024);
      m_scores.resize(m_index + 1024, HIGH_SCORE);
    }
    m_names[m_index] = stateName;
    m_scores[m_index] = score;
    m_index++;
  }


  virtual double GetScore(const string & stateName) const {
    int i;
    for (i=0; i<m_index; i++) {
      if (stateName == m_names[i]) {
	return m_scores[i];
      }
    }
    cout << "ERROR: No score for state " << stateName << endl;
    return HIGH_SCORE;
  }

 private:
  svec<string> m_names;
  svec<double> m_scores;
  int m_index;

};


//===========================================================
class HMMWord
{
 public:
  HMMWord() {
    HMMStateScore dummy;
    m_states.push_back(dummy);
    m_stateName.push_back("inscore");
    m_stay.push_back(0.);
    m_trans.push_back(0.);
  }

  void SetUp(const string & name) {
    m_name = name;
  }


  void SetStayPenalties(double d) {
    for (int i=1; i<m_stay.isize(); i++)
      m_stay[i] = d;
  }

  void AddState(const string & stateName, double stay = 0., double trans = 0.) {
    HMMStateScore dummy;
    m_states.push_back(dummy);
    m_stateName.push_back(stateName);
    m_stay.push_back(stay);
    m_trans.push_back(trans);
  }


  // Updating functions
  double GetOutScore() const {
    return m_states[m_states.isize()-1].GetScore() + m_trans[m_states.isize()-1];
  }
  int GetOutBack() const {
    return m_states[m_states.isize()-1].GetBack();
  }
  
  void SetIn(double score, int back) {
    m_states[0].Merge(score, back);
  }

  void Update(const HMMStateFrameScoreSource & a);


  const string & GetName() const {return m_name;}

 private:

  svec<HMMStateScore> m_states;
  svec<string> m_stateName;
  
  string m_name;
  svec<double> m_stay;
  svec<double> m_trans;
  
};



class HMMTraceback
{
 public:
  HMMTraceback() {}

  void GetBest(svec<HMMTracebackNode> & out,
	       const svec<HMMWord> & words,
	       const HMMHistory & h) const;

};





class HMMTransition
{
 public:
  HMMTransition() {
    m_default = 0.;
    m_bUseAcc = true;
  }

  void Resize(int n, double def = 100.) {
    m_n = n;
    m_matrix.resize(n * n, def);
    m_default = def;
    m_bUseAcc = true;
  }

  void Set(int from, int to, double val) {
    Get(from, to) = val;
    CheckDefault(val);
  }
 
  void SetAll(double d) {
    m_default = d;
    m_bUseAcc = true;
    for (int i=0; i<m_matrix.isize(); i++)
      m_matrix[i] = d;
  }

  double GetDefault() const {return m_default;}
  
  double & Get(int i, int j) {
    return m_matrix[m_n * i + j];
  }

  double GetTrans(int i, int j) const {
    return m_matrix[m_n * i + j];
  }

  bool IsAccelerated() const {return m_bUseAcc;}

  int GetSize() const {return m_n;}

 private:
  bool CheckDefault(double d) {
    double diff = d - m_default;
    if (diff < 0.01 && diff > -0.01)
      return true;
    m_bUseAcc = false;
    return false;
  }

  svec<double> m_matrix;
  double m_default;
  bool m_bUseAcc;

  int m_n;
};




//======================================================
class HMMUpdate
{
 public:

  HMMUpdate() {
    m_frame = 0;
  }


  int AddWord(const string & name) {
    HMMWord tmp;
    tmp.SetUp(name);
    m_words.push_back(tmp);
    m_trans.Resize(m_words.isize());
    return m_words.isize()-1;
  }

  void AddState(int wordIndex, const string & stateName, double stay = 0., double trans = 0.) {
    HMMWord & w = m_words[wordIndex];
    w.AddState(stateName, stay, trans);
  }

  void SetStayPenalty(double d) {
    for (int i=0; i<m_words.isize(); i++)
      m_words[i].SetStayPenalties(d);
  }

  void SetTransitionScore(int i, int j, double score) {
    m_trans.Set(i, j, score);
  }

  void SetAllTransitionScores(double score) {
    m_trans.SetAll(score);
  }

  void SetTransitionScore(const string & word1, const string & word2, double score);

  void Update(const HMMStateFrameScoreSource& a);


  void TraceBack(svec<HMMTracebackNode> & out) {
    HMMTraceback tb;

    tb.GetBest(out,
	       m_words,
	       m_history);

  }


  int GetWordCount() const {return m_words.isize();}

 private:
  int WordID(const string & s) {
    for (int i=0; i<m_words.isize(); i++) {
      if (m_words[i].GetName() == s)
	return i;
    }
    return -1;
  } 

  // The "vocabulary"
  svec<HMMWord> m_words;
  HMMHistory m_history;
  HMMTransition m_trans;
  int m_frame;
};




void Read(HMMUpdate & update, const string & config);





















#endif 



