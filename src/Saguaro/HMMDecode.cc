
#ifndef FORCE_DEBUG
#define NDEBUG
#endif

#include "src/Saguaro/HMMDecode.h"
#include "base/FileParser.h"


void Read(HMMUpdate & update, const string & config)
{

  int i;

  cout << "Setting up HMM..." << endl;
 
  FlatFileParser parserConfig;
  
  parserConfig.Open(config);

  while (parserConfig.ParseLine()) {
    if (parserConfig.GetItemCount() == 0)
      continue;
 
    const string & name = parserConfig.AsString(0);
 
    if (name == "all_transitions") {
      double score = parserConfig.AsFloat(1);
      update.SetAllTransitionScores(score);
      continue;
    }
    if (name == "stay_penalty") {
      double score = parserConfig.AsFloat(1);
      update.SetStayPenalty(score);
      continue;
    }

    if (name == "transition" || name == "Transition") {
      const string & word1 = parserConfig.AsString(1);
      const string & word2 = parserConfig.AsString(2);
      double score = parserConfig.AsFloat(3);
      update.SetTransitionScore(word1, word2, score);
      continue;
    }

    int n = update.AddWord(name);
    cout << "Adding word " << name << " as # " << n << endl;
    for (i=1; i<parserConfig.GetItemCount(); i++) {
      const string & state = parserConfig.AsString(i);
      update.AddState(n, state);
      cout << "     adding state " << state << endl;
    }
  }

}

void HMMWord::Update(const HMMStateFrameScoreSource & a) {
  int i;
  int n = m_stateName.isize();
  
  //cout << "Update word " << m_name << endl;

  for (i=n-1; i>0; i--) {    
    //cout << "State=" << i << endl;
    HMMStateScore & curr = m_states[i];
    HMMStateScore & prev = m_states[i-1];
    
    curr.AddScore(m_stay[i]);
    curr.Merge(prev.GetScore() + m_trans[i-1], prev.GetBack());      
    curr.AddScore(a.GetScore(m_stateName[i]));
    //cout << "Frame score=" << a.GetScore(m_stateName[i]) << endl;
  }
  //cout << "Updating word " << m_name << " score=" << m_states[1].GetScore() << endl;
  m_states[0].Reset();
}



//===============================================================

void HMMTraceback::GetBest(svec<HMMTracebackNode> & out,
			   const svec<HMMWord> & words,
			   const HMMHistory & h) const {
  int i;
  double low = h[h.GetLastFrameStart()].GetScore();
  int index = h.GetLastFrameStart();

  //cout << " start=" << h.GetLastFrameStart() << " end=" << h.GetLastFrameEnd() << endl;

  for (i=h.GetLastFrameStart(); i<h.GetLastFrameEnd(); i++) {
    //for (i=0; i<=h.GetLastFrameEnd(); i++) {
    //cout << i << "  score: " << h[i].GetScore() << " back=" << h[i].GetBack() << endl;
    if (h[i].GetScore() < low) {
      low = h[i].GetScore();
      index = i;
    }
  }
  
  // count them
  int n = 0;
  int index1 = index;
  while (index1 != -1) {      
    index1 = h[index1].GetBack();
    n++;
  }
  
  out.resize(n);
  
  
  index1 = index;
  while (index1 != -1) {      
    double score = h[index1].GetScore();
    int entry = h[index1].GetExit();
    const string & name = words[h[index1].GetWord()].GetName();
    
    index1 = h[index1].GetBack();
    n--;
    int exit = 0;
    if (index1 != -1)
      exit = h[index1].GetExit();
    
    out[n].Set(exit, entry, score, name, h[index1].GetWord());     
  }
  
}




void HMMUpdate::SetTransitionScore(const string & word1, const string & word2, double score)
{
  int k;
  if (word1 == "*") {
    int l = WordID(word2);
    cout << "Semi-explicit transition: " << word1 << " -> " << word2 << ": " << score << endl;
    for (k=0; k<m_trans.GetSize(); k++) {
      m_trans.Set(k, l, score);      
    }
    return;
  }
  if (word2 == "*") {
    int l = WordID(word1);
    cout << "Semi-explicit transition: " << word1 << " -> " << word2 << ": " << score << endl;
    for (k=0; k<m_trans.GetSize(); k++) {
      m_trans.Set(l, k, score);      
    }
  
    return;
  }
  int i = WordID(word1);
  int j = WordID(word2);
  cout << "Explicit transition: " << word1 << " -> " << word2 << ": " << score << endl;
  m_trans.Set(i, j, score);

}

//===========================================================
void HMMUpdate::Update(const HMMStateFrameScoreSource & a) {
  //m_history.PossiblySize(m_words.isize() * m_frames);
  int i, j;
  if (m_frame == 0) {
    for (i=0; i<m_words.isize(); i++) {
      m_words[i].SetIn(0., -1);
    }
  } else {

    if (false /*m_trans.IsAccelerated()*/) {
      //cout << "Running in accelerated mode." << endl;

      double bestOutScore = m_words[0].GetOutScore();
      int bestOutBack = -1;

      for (i=0; i<m_words.isize(); i++) {
	double outScore = m_words[i].GetOutScore();
	int outBack = m_words[i].GetOutBack();
	
	outBack = m_history.Add(m_frame-1, outBack, outScore, i);
	
	if (outScore < bestOutScore) {
	  bestOutBack = outBack;
	  bestOutScore = outScore;
	}

      }
      double inter = m_trans.GetDefault();
      for (j=0; j<m_words.isize(); j++) {     
	m_words[j].SetIn(bestOutScore + inter, bestOutBack);
      }
 

    } else {

      for (i=0; i<m_words.isize(); i++) {
	double outScore = m_words[i].GetOutScore();
	int outBack = m_words[i].GetOutBack();
	
	outBack = m_history.Add(m_frame-1, outBack, outScore, i);
	
	for (j=0; j<m_words.isize(); j++) {
	  double inter = m_trans.GetTrans(i, j);
	  m_words[j].SetIn(outScore + inter, outBack);
	}
      }
    }
    
  }
  // Update all words
  
  for (i=0; i<m_words.isize(); i++) {
    m_words[i].Update(a);
  }
  
  m_frame++;
}
