#include "src/Optical/OptiMap.h"
#define FORCE_DEBUG


double OpticalMap::Average()
{
  int i, j;
  double n = 0.;
  double d = 0.;
  for (i=0; i<isize(); i++) {
    OpticalSeq &s = m_data[i];
    for (j=1; j<s.isize();j++) {
      d += s[j]-s[j-1];
      n += 1.;
    }
  }
  return d/n;
}


double OptiAligner::AlignAuto(const OpticalSeq & a, const OpticalSeq & b)
{
  //cout << "ALIGN" << endl;
  double maxSite = a[a.isize()-1];
  if (b[b.isize()-1] > maxSite)
    maxSite = b[b.isize()-1];
 
  double scale = (m_size - 5)/maxSite;

  double auto_a = AlignInt(a, a, scale, 100000.);
  double auto_b = AlignInt(b, b, scale, 100000.);
 
  double thresh = auto_a;
  if (auto_b > auto_a)
    thresh = auto_b;
  thresh *= 0.5;

  
  double r = AlignInt(a, b, scale, thresh);
  //cout << "a=" << auto_a << " b=" << auto_b << " t=" << thresh << " r=" << r << endl;
  //cout << "Done." << endl;

  if (r > thresh)
    return r;
  else
    return 0;
}

double OptiAligner::Align(const OpticalSeq & a, const OpticalSeq & b, bool bPrint){
  double maxSite = a[a.isize()-1];
  if (b[b.isize()-1] > maxSite)
    maxSite = b[b.isize()-1];
 
  double scale = (m_size - 5)/maxSite;
  double thresh = 0.;
  if (!bPrint)
    thresh = 100000.;
  double r = AlignInt(a, b, scale, thresh);
  //double r = AlignAuto(a, b);
  return r;
}

double OptiAligner::AlignInt(const OpticalSeq & a, const OpticalSeq & b,
			   double scale, double thresh)
{
  svec<float> one, two;
  one.resize(m_size*2, 0.);
  two.resize(m_size*2, 0.);
  svec<float> out;
  out.resize(m_size*2, 0.);

  //cout << "Max=" << maxSite << " scale=" << scale << endl;; 
 

  Fit(one, a, scale);
  Fit(two, b, scale);

  m_xc.DoOne(out, one, two);

  int i;
  float max = 0;
  int maxPos = 0;
  for (i=0; i<out.isize(); i++) {
    //cout << i << " " << out[i] << endl;
    if (out[i] > max) {
      max = out[i];
      maxPos = i;
    }
  }
  //if (max < thresh)
  //return false;
  if (max > thresh)
    cout << a.Name() << " <-> " << b.Name() << " " << (double)(maxPos-m_size)/scale << "\t" << maxPos-m_size << "\t" << max << endl; 
  
  static int count = 0;
  if (max > thresh && a.Name() != b.Name()) {
    cout << "DETAILS:" << endl;
    
    int start = maxPos-m_size;
    if (start > 0) {
      for (i=0; i<m_size+start; i++) {
	if (i-start < 0 || i > m_size)
	  continue;
	cout << "ARRAY" << count << " "  << i << "\t" <<  two[i] << "\t" << i-start << "\t" << one[i-start] << endl;
      }
    } else {
      start = -start;
      for (i=0; i<m_size+start; i++) {
	if (i-start < 0 || i > m_size)
	  continue;
	cout << "ARRAY" << count << " " << i-start << "\t" <<  two[i-start] << "\t" << i << "\t" << one[i] << endl;
      }
    }
    count++;
  }
  return max;
}
 
double Func(double x, double x0, double a)
{
  double e = a*(x-x0);
  e *= e;
  double f = exp(-e);
  //cout << "e=" << e << " f=" << f << " x=" << x << " x0=" << x0 <<  endl;
  return f;
}

void OptiAligner::Fit(svec<float> & o, const OpticalSeq & a, double scale)
{
  int i, j;
  int win = 150;
  //cout << "Start" << endl;
  for (i=0; i<a.isize(); i++) {
    int p = (int)((double)a[i]*scale);
    //cout << "p=" << p << endl;
    if (p >= m_size) {
      cout << "ERROR" << endl;
      break;
    }
    //cout << "peak at " << p << endl;
    //double v = 1.;
    for (j=p-win; j<p+win; j++) {
      if (j >= m_size || j < 0)
	continue;
      double f = Func((double)j, p, m_decay);
      o[j] += f;
      if (o[j] > 1.)
	o[j] = 1.;
      //cout << "p=" << p << " j=" << j << " f=" << f << endl;
      //v *= m_decay;
    }
    /*
    v = m_decay;
    for (j=p-1; j>p-win; j--) {
      if (j < 0 || j >= m_size)
	break;
	o[j] += v;
      v *= m_decay;
      }*/
  }



  double sum = 0.;
  /* for (i=0; i<m_size; i++)
    sum += o[i];
  sum /= 64.;
  for (i=0; i<m_size; i++) {
    o[i] /= sum;
    }*/
 
  sum = 0;
  for (i=0; i<m_size; i++)
    sum += o[i];
  sum /= (double)m_size;
  for (i=0; i<m_size; i++) {
    o[i] -= sum;
    // cout << i << " " << o[i] << endl;
  }
  // cout << "End" << endl;

  //for (i=0; i<o.isize(); i++) {
  //  cout << i << " " << o[i] << endl;
  //}
}

