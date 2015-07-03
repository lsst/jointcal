// -*- C++ -*-
// 
// \file globalval.h
// 
// Last modified: $Date: 2010/09/09 18:10:42 $
// By:            $Author: seb $
// 

#ifndef GLOBALVAL__H
#define GLOBALVAL__H


#include <string>
#include <list>
#include <vector>
#include <map>
#include <iostream>


namespace lsst {
namespace meas {
namespace simastrom {



//! to store in files things like "Key value(s)" things.
// what exactly does this class bring except complicating a std::map?
// answer : use a std::map if you prefer.
class GlobalVal : private  std::map<std::string, std::vector<std::string> > {
public :

  GlobalVal() {};

  bool AddKey(const std::string &Key, const std::vector<std::string> &Values);
  
  bool AddKey(const std::string &Key, const std::string &Value);
  
  bool AddKey(const std::string& Key, const std::list<std::string>& Values);
  
  bool AddKey(const std::string &Key, const std::vector<double> &Values);
  
  bool AddKey(const std::string &Key, const double &Value);

  bool AddKey(const std::string& Key, const std::list<double>& Values);
  
  unsigned NKey() const;

  bool HasKey(const std::string &Key) const;
  
  std::string         getStringValue(const std::string& Key) const;
  
  std::vector<std::string> getStringValues(const std::string& Key) const;
  
  double         getDoubleValue(const std::string& Key) const;
  void           setDoubleValue(const std::string &Key, double val) ;
  std::vector<double> getDoubleValues(const std::string& Key) const;
  void           setDoubleValues(const std::string &Key, const std::vector<double>& vals);
  std::vector<std::string> OutputLines() const;

  void Erase(const std::string& Key);
  void Dump() const;
  void AppendTo(GlobalVal & glob) const ;
  void AppendTo(std::map<std::string, std::string>& globalKeys) const ;

  bool ProcessLine(const std::string &Line);

  // this constructor is to be used when you only need to read '@' lines in a file
  GlobalVal(const std::string &FileName);


private:
  template<typename T>
  bool GenericAddKey(const std::string& Key, const T& Values);
};

std::ostream& operator << (std::ostream &, const GlobalVal &G);

}}}

#endif /* GLOBALVAL__H */
