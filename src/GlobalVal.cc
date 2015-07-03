// 
// \file globalval.cc
// 
// Last modified: $Date: 2010/09/09 18:10:42 $
// By:            $Author: seb $
// 
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h> // for atof
#include <string.h>// for strlen
#include <cstdio>
#include <stdlib.h>

#include "lsst/meas/simastrom/GlobalVal.h"

namespace lsst {
namespace meas {
namespace simastrom {

using namespace std;

static int RemovePattern(string &Source, const string &aPattern)
{
  // cout << " Removing '" << aPattern << "' from '"<< Source << "'\n";
  int iter = 0;
  string::size_type pos = Source.find(aPattern);
  while (pos != Source.npos)
    {
      Source.erase(pos, aPattern.length()); 
      iter++;
      pos = Source.find(aPattern);
    }
  return iter;
}

static void DecomposeString(vector<string> &SubStrings, const string &Source, 
		     const char *tokens)
{
  string::size_type start = 0;
  string::size_type tokpos = Source.find_first_of(tokens,start);

  while (tokpos != Source.npos) 
    {
      string sub = Source.substr(start,tokpos-start);
      RemovePattern(sub," ");
      RemovePattern(sub,"\t");
      if (sub.length() > 0) SubStrings.push_back(sub);
      start = tokpos+1;
      tokpos = Source.find_first_of(tokens,start);
    }
  // get the last element
  string sublast = Source.substr(start, Source.length());
  RemovePattern(sublast," ");
  RemovePattern(sublast,"\t");
  if (sublast.length() > 0) SubStrings.push_back(sublast);

  // for (unsigned i=0; i<SubStrings.size(); ++i) 
  //   cout << i << "='"<< SubStrings[i] << "'"<<endl;
}



bool GlobalVal::HasKey(const string &Key) const
{
  return (find(Key) != end());
}


unsigned GlobalVal::NKey() const
{
  return this->size();
}


template<typename T>
bool GlobalVal::GenericAddKey(const string& Key, const T& Values) {
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
//      throw(PolokaException(" cannot have twice the same key val in GlobalVal "+Key));
      return false;
    }
    

  typename T::const_iterator I;
  for(I=Values.begin();I!=Values.end();I++) {
    stringstream sstrm;    
    sstrm << *I;
    (*this)[Key].push_back(sstrm.str());
  }
  return true;
}


bool GlobalVal::AddKey(const string& Key, const list<string>& Values)
{
  return GenericAddKey(Key, Values);
}


bool GlobalVal::AddKey(const string &Key, const vector<string> &Values)
{
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
//      throw(PolokaException(" cannot have twice the same key val in GlobalVal "+Key));
      return false;
    }
  
  (*this)[Key] = Values;
  return true;
}

bool GlobalVal::AddKey(const string &Key, const string &Value)
{
  vector<string> tmp;
  tmp.push_back(Value);
  return AddKey(Key,tmp);
}


bool GlobalVal::AddKey(const string &Key, const vector<double> &Values)
{
  if (HasKey(Key))
    {
      cerr << " cannot have twice the same key val in GlobalVal " << Key 
	   << endl;
//      throw(PolokaException(" cannot have twice the same key val in GlobalVal "+Key));
      return false;
    }
  
  unsigned int i;
  for(i=0;i<Values.size();i++) {
    stringstream sstrm;
    sstrm << Values[i];
    (*this)[Key].push_back(sstrm.str());
  }
  //  (*this)[Key] = Values;
  return true;
}


bool GlobalVal::AddKey(const string& Key, const list<double>& Values)
{
  return GenericAddKey(Key, Values);
}


bool GlobalVal::AddKey(const string &Key, const double &Value)
{
  vector<double> tmp;
  tmp.push_back(Value);
  return AddKey(Key,tmp);
}

void GlobalVal::Erase(const string &Key)
{
  map<string, vector<string> >::iterator it = find(Key);
  if(it!= end())
    erase(it);
}


string GlobalVal::getStringValue(const string& Key) const
{
  const_iterator i = find(Key);
  if (i == end())
    {
      cerr << " could not find key " << Key << endl;
//      throw(PolokaException(" could not find key "+Key));
      return "";
    }
  else return i->second[0];
}


vector<string> GlobalVal::getStringValues(const string& Key) const
{
  const_iterator i = find(Key);
  if(i == end())
    {
      cerr << " could not find key " << Key << endl;
//      throw(PolokaException(" could not find key "+Key));
    }
  return i->second;
}


double GlobalVal::getDoubleValue(const string &Key) const
{
  const_iterator i = find(Key);
  if (i == end())
    {
      cerr << " could not find key " << Key << endl;
//      throw(PolokaException(" could not find key "+Key));
    }
  else return atof(i->second[0].c_str());
}

void GlobalVal::setDoubleValue(const string &Key, double val) 
{
  vector<string>& vals=(*this)[Key];
  vals.clear();
  char c[100] ;
  sprintf(c,"%f",val);
  vals.push_back(c);
}

void GlobalVal::setDoubleValues(const string &Key, const vector<double>& Vals) 
{
  vector<string>& vals=(*this)[Key];
  vals.clear();
  for (size_t iv=0; iv<Vals.size(); iv++) {
    char c[100] ;
    sprintf(c,"%f",Vals[iv]);
    vals.push_back(c);
  }
}


vector<double> GlobalVal::getDoubleValues(const string &Key) const
{
  const_iterator i = find(Key);
  if (i == end())
    {
      cerr << " could not find key " << Key << endl;
      //      throw(PolokaException(" could not find key "+Key));
    }
  vector<double> ret;
  for(unsigned int k=0;k<i->second.size();k++)
    ret.push_back(atof(i->second[k].c_str()));
  return ret;
}

void GlobalVal::AppendTo(GlobalVal & glob) const
{
for (const_iterator i = begin(); i != end(); ++i)
    {
      //      const vector<double> &values = i->second;
      const vector<string> &values = i->second;
      if (values.size() == 0) glob.AddKey( i->first, "") ;
      else
	glob.AddKey( i->first, values) ;
    }
}


void GlobalVal::AppendTo(map<string, string>& globalKeys) const
{
for (const_iterator i = begin(); i != end(); ++i)
    {
      //      const vector<double> &values = i->second;
      const vector<string> &values = i->second;
      if (values.size() == 0) globalKeys[i->first]=""; 
      else
	{
	  string l_values ;
	  for(int ii = 0 ; ii < values.size(); ii++)
	    l_values += " " + values[ii];
	  globalKeys[i->first]=l_values ;
	}
    }
}






void GlobalVal::Dump() const
{
  for (const_iterator i = begin(); i != end(); ++i)
    {
      const vector<string> &values = i->second;
      //if (values.size() == 0) continue;
      cerr << i->first;
      for (unsigned k=0; k < values.size(); ++k) cerr << ' ' << values[k];
      cerr << endl ;
    }
  return ;
}

vector<string> GlobalVal::OutputLines() const
{
  vector<string> out;
  for (const_iterator i = begin(); i != end(); ++i)
    {
      //      const vector<double> &values = i->second;
      const vector<string> &values = i->second;
      if (values.size() == 0) continue;
      ostringstream s;
      s << i->first;
      for (unsigned k=0; k < values.size(); ++k) s << ' ' << values[k];
      out.push_back(s.str());
    }
  return out;
}



GlobalVal::GlobalVal(const std::string &FileName)
{
  std::ifstream r(FileName.c_str());
  char c ;
  char buff[4096];
  while( r >> c ) // to test eof
    {
      r.unget() ;
      if ( (c == '@') ) 
	{
	  r.getline(buff,4096); 
	  ProcessLine(buff);
	  continue;
	}
      if ( (c=='#') || isdigit(c)) break;
    }
  r.close();
}




#ifdef STORAGE
static void explode(vector<string>& v, string const& str, char tok=',')
{
  string s, sc=str;
  string::size_type /* p0=0, */ p1;

  if(str == "*") {
    stringstream sstrm;
    for(int i=0;i<5;i++) {
      sstrm << i;
      v.push_back(sstrm.str());
      sstrm.str("");
    }
    return;
  }

  do {
    p1 = sc.find_first_of(tok);
    s = sc.substr(0, p1);
    v.push_back(s);
    if(p1 != string::npos)
      sc = sc.substr(p1+1, string::npos);
  } while(p1 != string::npos);

}
#endif


bool GlobalVal::ProcessLine(const string &Line)
{
  // use standard C stuff because it enables error checking
  const char *s = Line.c_str();
  if (*s == '@') s++;

  char buf[256];
  if (sscanf(s,"%s",buf)!= 1)
    {
      cerr << " could no read a key in " << Line << endl;
      //      throw(PolokaException(" could no read a key in "+Line));
      return false;
    }
  s += strlen(buf);
  string key(buf);
  if (HasKey(key))
    {
      cerr <<" already have a key labelled " << key << endl;
      //      throw(PolokaException(" already have a key labelled "+key+" in line :\n"+Line));
    }
  vector<string> &values = (*this)[key];
  string strtmp = s;
  DecomposeString(values, strtmp, " ");
  return true;
}

ostream& operator << (ostream &S , const GlobalVal &G)
{
  vector<string> output(G.OutputLines());
  for (unsigned k=0; k< output.size(); ++k) S << output[k] << endl;
  return S;
}

}}}
