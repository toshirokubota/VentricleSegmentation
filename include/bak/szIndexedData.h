#ifndef ___INDEXED_DATA_H____
#define ___INDEXED_DATA_H____

class indexedData  
{
public:
  double data;
  int index;
};

inline bool 
operator<(const indexedData& s1, const indexedData& s2) {
  return s1.data < s2.data;
}

#endif /* ___INDEXED_DATA_H____ */