#ifndef _GEN_KDTREE_H_
#define _GEN_KDTREE_H_

namespace PG
{

template<int N>
struct KDBox {
 double bmin[N], bmax[N];
 int id;
};

template<int N>
struct KDTNode {
 KDTNode* l;
 KDTNode* r;
 std::vector<KDBox<N>*> bs;
 double bmin[N], bmax[N];
 inline static KDTNode* Build(std::vector<KDBox<N>*>& bs);
 inline static bool Hit(
  KDTNode* node, const double* p0, const double* p1, std::vector<int>& ids);
 inline KDTNode();
 inline ~KDTNode();
};

template<int N>
KDTNode<N>::KDTNode() : l(nullptr), r(nullptr) {}

template<int N>
KDTNode<N>::~KDTNode()
{
 if (l) delete l;
 if (r) delete r;
}

template<int N>
KDTNode<N>* KDTNode<N>::Build(std::vector<KDBox<N>*>& boxes)
{
 KDTNode* node = new KDTNode();
 for (int i = 0; i < N; ++i) node->bmin[i] = 0.0;
 for (int i = 0; i < N; ++i) node->bmax[i] = 0.0;
 node->bs = boxes;

 if (boxes.empty()) return node;

 for (int i = 0; i < N; ++i) node->bmin[i] = boxes[0]->bmin[i];
 for (int i = 0; i < N; ++i) node->bmax[i] = boxes[0]->bmax[i];

 if (boxes.size() == 1)
 {
  node->l = new KDTNode();
  node->r = new KDTNode();
  node->l->bs = std::vector<KDBox<N>*>();
  node->r->bs = std::vector<KDBox<N>*>();
  return node;
 }

 const int ntris = (int) boxes.size();
 for (int i = 1; i < ntris; ++i)
 {
  for (int j = 0; j < N; ++j)
   if (node->bmin[j] > boxes[i]->bmin[j])
    node->bmin[j] = boxes[i]->bmin[j];
  for (int j = 0; j < N; ++j)
   if (node->bmax[j] < boxes[i]->bmax[j])
    node->bmax[j] = boxes[i]->bmax[j];
 }

 double mid[N]; for (int i = 0; i < N; ++i) mid[i] = 0;
 for (int i = 0; i < ntris; ++i)
  for (int j = 0; j < N; ++j)
   mid[j] += (boxes[i]->bmin[j] + boxes[i]->bmax[j]) * 0.5;
 for (int i = 0; i < N; ++i)
  mid[i] /= (double) ntris;

 int axis = 0;
 for (int j = 1; j < N; ++j)
  if ((node->bmax[axis] - node->bmin[axis]) <
      (node->bmax[j]    - node->bmin[j]) )
   axis = j;

 std::vector<KDBox<N>*> left, right;
 for (int i = 0; i < ntris; ++i)
  (mid[axis] >= (boxes[i]->bmin[axis] + boxes[i]->bmax[axis]) * 0.5) ?
   right.push_back(boxes[i]) : left.push_back(boxes[i]);

 if (left.empty() && !right.empty()) left = right;
 if (right.empty() && !left.empty()) right = left;

 int matches(0);
 const int nleft = (int) left.size();
 const int nright = (int) right.size();
 for (int i = 0; i < nleft; ++i)
  for (int j = 0; j < nright; ++j)
   if (left[i] == right[j])
    ++matches;

 if (((matches / (float) nleft) < 0.5f) &&
     ((matches / (float) nright) < 0.5f))
 {
  node->l = Build(left);
  node->r = Build(right);
 }
 else
 {
  node->l = new KDTNode<N>();
  node->r = new KDTNode<N>();
  node->l->bs = std::vector<KDBox<N>*>();
  node->r->bs = std::vector<KDBox<N>*>();
 }

 return node;
}

template<int N>
static bool AABBRay(
 const double* bmin, const double* bmax, const double* p0, const double* p1)
{
 double c[N]; for (int i=0; i<N; ++i) c[i] = 0;
 for (int i = 0; i < N; ++i) c[i] = (bmin[i] + bmax[i]) * 0.5; // Box center-point
 double e[N]; for (int i=0; i<N; ++i) e[i] = 0;
 for (int i = 0; i < N; ++i) e[i] = bmax[i] - c[i]; // Box halflength extents
 double m[N]; for (int i=0; i<N; ++i) m[i] = 0;
 for (int i = 0; i < N; ++i) m[i] = (p0[i] + p1[i]) * 0.5; // Segment midpoint
 double d[N]; for (int i=0; i<N; ++i) d[i] = 0;
 for (int i = 0; i < N; ++i) d[i] = p1[i] - m[i]; // Segment halflength vector
 for (int i = 0; i < N; ++i) m[i] = m[i] - c[i]; // Translate box and segment to origin
 // Try world coordinate axes as separating axes
 double ad[N];
 for (int i = 0; i < N; ++i) ad[i] = fabs(d[i]);
 for (int i = 0; i < N; ++i)
  if (fabs(m[i]) > e[i] + ad[i])
   return false;
 // Add in an epsilon term to counteract arithmetic errors when segment is
 // (near) parallel to a coordinate axis
 for (int i = 0; i < N; ++i)
  ad[i] += std::numeric_limits<double>::epsilon();
 // Try cross products of segment direction vector with coordinate axes
 if (fabs(m[1] * d[2] - m[2] * d[1]) > e[1] * ad[2] + e[2] * ad[1]) return false;
 if (fabs(m[2] * d[0] - m[0] * d[2]) > e[0] * ad[2] + e[2] * ad[0]) return false;
 if (fabs(m[0] * d[1] - m[1] * d[0]) > e[0] * ad[1] + e[1] * ad[0]) return false;
 // No separating axis found; segment must be overlapping AABB
 return true;
}

template<int N>
bool KDTNode<N>::Hit(
 KDTNode* node, const double* p0, const double* p1, std::vector<int>& ids)
{
 if (AABBRay<N>(node->bmin, node->bmax, p0, p1))
 {
  if (!node->l->bs.empty() || !node->r->bs.empty())
  {
   bool hleft = Hit(node->l, p0, p1, ids);
   bool hright = Hit(node->r, p0, p1, ids);
   return hleft || hright;
  }
  else // leaf
  {
   bool hit(false);
   const int ntris = (int) node->bs.size();
   for (int i = 0; i < ntris; ++i)
    if (AABBRay<N>(node->bs[i]->bmin, node->bs[i]->bmax, p0, p1))
    {
     hit = true;
     ids.push_back(node->bs[i]->id);
    }
   return hit;
  }
 }
 return false;
}

}

#endif /* _GEN_KDTREE_H_ */