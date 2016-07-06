#ifndef ___SZ_GEOMETRY_H___
#define ___SZ_GEOMETRY_H___
#include <iostream>
#include <list>
#include <vector>
using namespace std;

class szEdge; //forward declaration
class szFace; //forward declaration

class szVertex {
public:
  szVertex(long x=0, long y=0, long z=0) {
    //, szEdge* pedge=NULL) {
    m_lX = x;
    m_lY = y;
    m_lZ = z;
    ///m_peIncident = pedge;
  }

  ~szVertex() {
    //m_peIncident = NULL;
  }

  long getX() const {return m_lX;}
  long getY() const {return m_lY;}
  long getZ() const {return m_lZ;}
  //szEdge* getIncidentEdge() const {return m_peIncident;}
  void setX(long x) {m_lX=x;}
  void setY(long y) {m_lY=y;}
  void setZ(long z) {m_lZ=z;}
  //void setIncidentEdge(szEdge* pedge) {m_peIncident=pedge;}

protected:
  long m_lX;
  long m_lY;
  long m_lZ;
  //szEdge* m_peIncident;
};

class szEdge {
public:
  szEdge(szVertex* porigin=NULL, 
    szEdge* pe_twin=NULL, szEdge* pe_next=NULL, 
    szEdge* pe_prev=NULL, szFace* pf_inc=NULL) {
    m_vOrigin = porigin;
    m_eTwin = pe_twin;
    m_eNext = pe_next;
    m_ePrev = pe_prev;
    m_fFace = pf_inc;
  }

  ~szEdge() {
    m_vOrigin = NULL;
    m_eTwin = NULL;
    m_fFace = NULL;
    m_eNext = NULL;
    m_ePrev = NULL;
  }

  szVertex* getOrigin() const {return m_vOrigin;}
  szEdge* getTwin() const {return m_eTwin;}
  szEdge* getNext() const {return m_eNext;}
  szEdge* getPrev() const {return m_ePrev;}
  szFace* getIncidentFace() const {return m_fFace;}
  void setOrigin(szVertex* pv) {m_vOrigin=pv;}
  void setTwin(szEdge* pedge) {m_eTwin=pedge;}
  void setNext(szEdge* pedge) {m_eNext=pedge;}
  void setPrev(szEdge* pedge) {m_ePrev=pedge;}
  void setIncidentFace(szFace* pface) {m_fFace=pface;}  

protected:
  szVertex* m_vOrigin;
  szEdge* m_eTwin;
  szFace* m_fFace;
  szEdge* m_eNext;
  szEdge* m_ePrev;
};


class szFace {
public:
  szFace() {
    m_peOuterComponent = NULL;
    m_peInnerComponent = NULL;
    //printf("Face %x allocated\n", (int)this);
  }

  /*szFace(szEdge* pedge=NULL, bool outer=true) {
    if(outer) {
      m_peOuterComponent = pedge;
      m_peInnerComponent = NULL;
    }
    else {
      m_peOuterComponent = NULL;
      m_peInnerComponent = pedge;
    }
  }*/

  ~szFace() {
    //if((int)this == 0x00520680)
    //  printf("Face %x deleted\n", (int)this);
    szEdge* peFirst=NULL;
    if(m_peOuterComponent!=NULL) 
      peFirst = m_peOuterComponent;
    else if(m_peInnerComponent != NULL) 
      peFirst = m_peInnerComponent;
    //else {
      //it has been manually cleaned up
    //  return;
    //}

    szEdge* pe = peFirst;
    szEdge* peNext;
    do {
      peNext = pe->getNext();
      delete pe;
      pe = peNext;
    }
    while(pe != peFirst); // && pe != NULL);
    m_peOuterComponent = NULL;
    m_peInnerComponent = NULL;
  }

  /*szFace(const szFace& face) {
    m_peOuterComponent = face.getOuterComponent();
    m_peInnerComponent = face.getInnerComponent();
  }

  const szFace& operator =(const szFace& face) {
    m_peOuterComponent = face.getOuterComponent();
    m_peInnerComponent = face.getInnerComponent();
    return *this;
  }*/

  szEdge* getOuterComponent() const {return m_peOuterComponent;}
  szEdge* getInnerComponent() const {return m_peInnerComponent;}
  void setOuterComponent(szEdge* pedge) {m_peOuterComponent = pedge;}
  void setInnerComponent(szEdge* pedge) {m_peInnerComponent = pedge;}

  vector<int> getNormal() const {
    szEdge* pe = m_peOuterComponent;
    szVertex* pv = pe->getOrigin();
    int x1=pv->getX();
    int y1=pv->getY();
    int z1=pv->getZ();
    pe = pe->getNext();
    pv = pe->getOrigin();
    int x2=pv->getX()-x1;
    int y2=pv->getY()-y1;
    int z2=pv->getZ()-z1;
    pe = pe->getNext();
    pv = pe->getOrigin();
    int x3=pv->getX()-x1;
    int y3=pv->getY()-y1;
    int z3=pv->getZ()-z1;
    
    vector<int> vnorm;
    vnorm.push_back(y2*z3-z2*y3);
    vnorm.push_back(z2*x3-x2*z3);
    vnorm.push_back(x2*y3-y2*x3);

    return vnorm;
  }

protected:
  szEdge* m_peOuterComponent;
  szEdge* m_peInnerComponent;
};

class szConflictGraph {
public:
  szConflictGraph() {
  }

  ~szConflictGraph() {
    m_lpVertex.clear();
    m_lpFace.clear();
  }

  list<szVertex*> m_lpVertex;
  list<szFace*> m_lpFace;
};


//double
//getArea(szVertex* p1, szVertex* p2, szVertex* p3);

double
getVolume(const szVertex* p1, const szVertex* p2, const szVertex* p3, const szVertex* p4);

bool
IsVisible(const szVertex* pv, const szFace* pf);

//bool
//IsCCW(szVertex* p1, szVertex* p2, szVertex* p3);

bool
CoLinear(const szVertex* p1, const szVertex* p2, const szVertex* p3);

bool
CoLinear(const szEdge* pe1, const szEdge* pe2);

bool
CoPlaner(const szVertex* p1, const szVertex* p2, const szVertex* p3, const szVertex* p4);

bool
CoPlaner(const szVertex* p1, const szFace* pf);

szFace*
MakeNewFace(szVertex* pv,
            szEdge* pedge);

szFace*
MakeNewFace(szVertex* pv1,
            szVertex* pv2,
            szVertex* pv3);

bool
IsTwin(const szEdge* pe1, const szEdge* pe2);

szEdge*
FindTwin(const szFace* pf, const szEdge* pe);

bool
MergeFaces(szFace* pf1, const szFace* pf2);

bool
RemoveColinearEdges(szFace* pf);

ostream& operator<<(ostream& s, const szVertex& v);
ostream& operator<<(ostream& s, const szEdge& v);
ostream& operator<<(ostream& s, const szFace& v);
ostream& operator<<(ostream& s, const szConflictGraph& v);
void PrintFace(const szFace& face, ostream& out);

#endif /* ___SZ_GEOMETRY_H___ */