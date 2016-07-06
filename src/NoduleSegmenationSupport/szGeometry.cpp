#include <math.h>
#include <szGeometry.h>

double
getArea(const szVertex* p1, 
        const szVertex* p2, 
        const szVertex* p3) 
{
  int dx1=p2->getX()-p1->getX();
  int dy1=p2->getY()-p1->getY();
  int dz1=p2->getZ()-p1->getZ();
  int dx2=p3->getX()-p1->getX();
  int dy2=p3->getY()-p1->getY();
  int dz2=p3->getZ()-p1->getZ();
  
  int ovx = dy1*dz2-dz1*dy2;
  int ovy = dz1*dx2-dx1*dz2;
  int ovz = dx1*dy2-dy1*dx2;  

  double area = sqrt((double)ovx*ovx+ovy*ovy+ovz*ovz)/2.0;
  return area;
}

double
getVolume(const szVertex* p1, 
          const szVertex* p2, 
          const szVertex* p3, 
          const szVertex* p4) 
{
  int dx1=p2->getX()-p1->getX();
  int dy1=p2->getY()-p1->getY();
  int dz1=p2->getZ()-p1->getZ();
  int dx2=p3->getX()-p1->getX();
  int dy2=p3->getY()-p1->getY();
  int dz2=p3->getZ()-p1->getZ();
  
  int ovx = dy1*dz2-dz1*dy2;
  int ovy = dz1*dx2-dx1*dz2;
  int ovz = dx1*dy2-dy1*dx2;  
  int dx3=p4->getX()-p1->getX();
  int dy3=p4->getY()-p1->getY();
  int dz3=p4->getZ()-p1->getZ();

  double vol = ((double)ovx*dx3+ovy*dy3+ovz*dz3)/6.0;
  return vol;
}

bool
IsVisible(const szVertex* pv, const szFace* pf) 
{
  vector<int> vn = pf->getNormal();
  szVertex* po = pf->getOuterComponent()->getOrigin();
  int x=pv->getX()-po->getX();
  int y=pv->getY()-po->getY();
  int z=pv->getZ()-po->getZ();
  int vol = x*vn[0]+y*vn[1]+z*vn[2];
  if(vol>0)
    return true;
  else
    return false;
}

bool
CoLinear(const szVertex* p1, 
         const szVertex* p2, 
         const szVertex* p3) 
{
  if(getArea(p1,p2,p3)==.0)
    return true;
  else
    return false;
}

bool
CoLinear(const szEdge* pe1, 
         const szEdge* pe2) 
{
  if(getArea(pe1->getOrigin(),pe1->getNext()->getOrigin(),pe2->getOrigin())==.0) 
  {
    if(getArea(pe1->getOrigin(),pe1->getNext()->getOrigin(),pe2->getNext()->getOrigin())==.0)
      return true;
  }
  return false;
}

bool
CoPlaner(const szVertex* p1, 
         const szVertex* p2, 
         const szVertex* p3, 
         const szVertex* p4) 
{
  if(getVolume(p1,p2,p3,p4)==.0)
    return true;
  else
    return false;
}

/*
wrapper for CoPlaner with a face and a vertex as arguements.
It simply gets three vertices from the face and call ColPlaner.
*/
bool
CoPlaner(const szVertex* p1, 
         const szFace* pf) 
{
  szEdge* pe = pf->getOuterComponent();
  if(pe==NULL) {
    pe = pf->getInnerComponent();
    if(pe == NULL)
      return false;
  }

  szVertex* pv[3];
  for(int i=0; pe!=NULL && i<3; ++i) {
    pv[i] = pe->getOrigin();
    pe = pe->getNext();
  }

  return CoPlaner(p1,pv[0],pv[1],pv[2]);
}

/*
Make a face from one vertex and one edge.
We do not have to check the correct order, as it should have been done when the edge 
pointed by pedge was created.
*/
szFace*
MakeNewFace(szVertex* pv,
            szEdge* pedge) 
{
  //there are three new edges
  szEdge* pedge1 = new szEdge;
  szEdge* pedge2 = new szEdge;
  szEdge* pedge3 = new szEdge;
  szFace* pface = new szFace;

  pedge1->setOrigin(pedge->getOrigin());
  pedge1->setNext(pedge2);
  pedge1->setPrev(pedge3);
  pedge1->setTwin(pedge->getTwin());
  pedge1->setIncidentFace(pface);

  pedge2->setOrigin(pedge->getNext()->getOrigin());
  pedge2->setNext(pedge3);
  pedge2->setPrev(pedge1);
  pedge2->setIncidentFace(pface);

  pedge3->setOrigin(pv);
  pedge3->setNext(pedge1);
  pedge3->setPrev(pedge2);
  pedge3->setIncidentFace(pface);

  pface->setOuterComponent(pedge1);

  return pface;
}

/*
Make a face from three vertices
Vertices need to be given in counter-clock wise when the surface is viewed from the outside
*/
szFace*
MakeNewFace(szVertex* pv1,
            szVertex* pv2,
            szVertex* pv3) 
{
  //there are six new edges
  szEdge* pedge1 = new szEdge;
  szEdge* pedge2 = new szEdge;
  szEdge* pedge3 = new szEdge;
  szFace* pface = new szFace;

  pedge1->setOrigin(pv1);
  pedge1->setNext(pedge2);
  pedge1->setPrev(pedge3);
  pedge1->setIncidentFace(pface);

  pedge2->setOrigin(pv2);
  pedge2->setNext(pedge3);
  pedge2->setPrev(pedge1);
  pedge2->setIncidentFace(pface);

  pedge3->setOrigin(pv3);
  pedge3->setNext(pedge1);
  pedge3->setPrev(pedge2);
  pedge3->setIncidentFace(pface);

  pface->setOuterComponent(pedge1);

  return pface;
}

/*
find the twin edge of a given edge in a given face.
*/
bool
IsTwin(const szEdge* pe1, 
       const szEdge* pe2) 
{
  if(pe1->getOrigin()==pe2->getNext()->getOrigin() &&
    pe2->getOrigin()==pe1->getNext()->getOrigin()) { //NOTE: the vertex order is reversed
    return true;
  }
  else
    return false;
}

/*
find the twin edge of a given edge in a given face.
*/
szEdge*
FindTwin(const szFace* pf, 
         const szEdge* pe) 
{
  szVertex* pv1 = pe->getOrigin();
  szVertex* pv2 = pe->getNext()->getOrigin();

  szEdge* pe2 = pf->getOuterComponent();
  if(pe2==NULL) {
    pe2 = pf->getInnerComponent();
    if(pe2==NULL)
      return NULL;
  }

  szEdge* peBegin = pe2;
  do {
    if(pe2->getOrigin()==pv2 && pe2->getNext()->getOrigin()==pv1) { //NOTE: the vertex order is reversed
      return pe2;
    }
    pe2 = pe2->getNext();
  }
  while(pe2!=peBegin && pe2!=NULL);

  return NULL;
}

//merge two faces (pf1 consumes pf2) along a pair of edges that are twin to each other.
//we assume that there is only one section that make contact with each other.  The section
//can be comprised of multiple edges, however.
bool
MergeFaces(szFace* pf1, 
           const szFace* pf2) 
{
  szEdge* peTwin1=pf1->getOuterComponent();
  szEdge* peBegin=peTwin1;
  szEdge* peTwin2=NULL;

  //find a twin pair
  do {
    peTwin2=FindTwin(pf2,peTwin1);
    if(peTwin2==NULL)
      peTwin1=peTwin1->getNext();
  }
  while(peTwin2==NULL && peTwin1!=peBegin);

  if(peTwin2==NULL)
    return false;

  //check how long the two surface touches
  szEdge* peMergeEnd1=peTwin1;
  szEdge* peMergeBegin2=peTwin2;
  while(IsTwin(peMergeEnd1->getNext(),peMergeBegin2->getPrev())) {
    peMergeEnd1=peMergeEnd1->getNext();
    peMergeBegin2=peMergeBegin2->getPrev();
  }
  szEdge* peMergeBegin1=peTwin1;
  szEdge* peMergeEnd2=peTwin2;
  while(IsTwin(peMergeBegin1->getPrev(),peMergeEnd2->getNext())) {
    peMergeBegin1=peMergeBegin1->getPrev();
    peMergeEnd2=peMergeEnd2->getNext();
  }

  //go around pf2 and copy edges
  szEdge* pe1=peMergeBegin1->getPrev();
  szEdge* pe2=peMergeEnd2->getNext();
  while(pe2!=peMergeBegin2) {
    szEdge* peNew = new szEdge;
    peNew->setOrigin(pe2->getOrigin());
    peNew->setIncidentFace(pf1);
    peNew->setPrev(pe1);
    peNew->setTwin(pe2->getTwin());
    pe1->setNext(peNew);
    pe1=peNew;
    pe2=pe2->getNext();
  }
  //make the final bridge between the two faces.
  pe1->setNext(peMergeEnd1->getNext());
  pe1->getNext()->setPrev(pe1);

  szEdge* peDel=peMergeBegin1;
  do {
    szEdge* pe=peDel->getNext();
    if(peDel==pf1->getOuterComponent())
      pf1->setOuterComponent(pe1);
    delete peDel;
    peDel=pe;
  }
  while(peDel!=pe1->getNext());
  
  return true;
}

bool
RemoveColinearEdges(szFace* pf) 
{
  szFace* gFace=(szFace*)0x4d1c90;

  szEdge* pe = pf->getOuterComponent();
  szEdge* peBegin=pe;
  bool first=true; //this flag takes care of the case when the first edge is co-linear
                   //without this, the do-while loop will not visit all the edges.
  do {
    if(CoLinear(pe->getOrigin(),pe->getNext()->getOrigin(),pe->getNext()->getNext()->getOrigin())) {
      szEdge* peDel = pe->getNext();
      pe->setNext(pe->getNext()->getNext());
      pe->getNext()->setPrev(pe);
      if(pe->getTwin()==0 && peDel->getTwin()!=NULL)
        pe->setTwin(peDel->getTwin());
      if(peDel->getTwin()!=NULL)
        peDel->getTwin()->setTwin(pe);
      delete peDel;
      if(peDel==peBegin) {//we just deleted the first one
        pf->setOuterComponent(pe);
        break; //important - the exit condition no longer valid
      }
    }
    else {
      first = false;
      pe=pe->getNext();
    }
  }
  while(first || pe!=peBegin);

  return true;
}


ostream& operator<<(ostream& s, const szVertex& v) {
  s << "(" << v.getX() << ", " << v.getY() << ", " << v.getZ() << ")";
  return s;
}

ostream& operator<<(ostream& s, const szFace& f) {
  szEdge* pe=f.getOuterComponent();
  szEdge* peBegin=pe;

  s << "F::";
  while(pe!=NULL) {
    s << *(pe->getOrigin());
    pe=pe->getNext();
    if(pe==peBegin)
      break;
    else
      s << " - ";
  }
  //s << endl;

  return s;
}

ostream& operator<<(ostream& s, const szEdge& f) {

  s << "E::";
  s << " " << *(f.getOrigin());
  s << " - " << *(f.getNext()->getOrigin());
  //s << " [" << f.getPrev();
  //s << " - " << f.getNext();
  //s << " + " << f.getTwin();
  //s << " = " << f.getIncidentFace();
  //s << "]" << endl;*/

  return s;
}

ostream& operator<<(ostream& s, const szConflictGraph& f) {
  list<szVertex*>::const_iterator lv;
  list<szFace*>::const_iterator lf;
  int count=0;
  for(lv=f.m_lpVertex.begin(), lf=f.m_lpFace.begin(); lv!=f.m_lpVertex.end(); ++lv, lf) {
    s << ++count << ": " << **lv << endl;
    s << "\t" << **lf << endl;
  }

  return s;
}

void
PrintFace(const szFace& face, ostream& out) {
  szEdge* pe=face.getOuterComponent();
  if(pe==NULL)
    return;

  out << "PrintFace:\n";
  szEdge* peBegin=pe;
  do {
    out << *(pe->getOrigin()) << " - ";

    pe=pe->getNext();
  }
  while(pe!=peBegin);
  out << endl;
}
