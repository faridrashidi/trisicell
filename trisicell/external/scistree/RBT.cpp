#include "RBT.h"

//////////////////////////////////////////////////////////////////////////////

// useful stuff
int GetNumRBT(int nlv) {
  int res = 1;
  for (int nr = 2; nr < nlv; ++nr) {
    //
    res *= 2 * nr - 1;
  }
  return res;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

int RBTNode ::idNodeNextToUse = 20000;

RBTNode ::RBTNode(RBTNode *pLeftParam, RBTNode *pRightParam)
    : pLeft(pLeftParam), pRight(pRightParam), pParent(NULL) {
  YW_ASSERT_INFO(pLeft != NULL && pRight != NULL, "Can not be NULL");

  // ensure children's parent are set
  pLeft->SetParent(this);
  pRight->SetParent(this);
  lvid = idNodeNextToUse++;
  SetHeight(-1.0);
}

// operation
RBTNode *RBTNode ::CopySubTree() {
  // copy entire subtree under it
  if (IsLeaf() == false) {
    // copy left/right subtrees
    RBTNode *pLT = pLeft->CopySubTree();
    RBTNode *pRT = pRight->CopySubTree();
    RBTNode *pNewNode = new RBTNode(pLT, pRT);
    // cout << "copy a internal node " <<", newnode = " << (int) pNewNode <<
    // endl;
    return pNewNode;
  } else {
    // copy self only
    RBTNode *pNewNode = new RBTNode(this->lvid);
    pNewNode->SetHeight(-1.0);
    // cout << "copy a leaf node "<< ", lvid = " << lvid  <<", newnode = " <<
    // (int) pNewNode << endl;
    return pNewNode;
  }
}

void RBTNode ::AddToLeftEdge(int lvidParam) {
  // ensure this is not a leaf
  YW_ASSERT_INFO(IsLeaf() == false, "Can not be a leaf");

  RBTNode *pInternal = pLeft->AddSibling(lvidParam);
  pInternal->SetParent(this);
  this->SetLeftChild(pInternal);
}

void RBTNode ::AddToRightEdge(int lvidParam) {
  YW_ASSERT_INFO(IsLeaf() == false, "Can not be a leaf");
  RBTNode *pInternal = pRight->AddSibling(lvidParam);
  pInternal->SetParent(this);
  this->SetRightChild(pInternal);
}

RBTNode *RBTNode ::AddSibling(int lvidParam) {
  // create a new sibling and a root
  RBTNode *pOther = new RBTNode(lvidParam);
  // cout << "Adding a leaf " << (int) pOther << endl;
  int existId = GetMinLeaveId();
  if (existId < lvidParam) {
    //
    RBTNode *pParent = new RBTNode(this, pOther);
    // cout << "Adding a node " << (int) pParent << endl;
    return pParent;
  } else {
    RBTNode *pParent = new RBTNode(pOther, this);
    // cout << "Adding a node " << (int) pParent << endl;
    return pParent;
  }
}

void RBTNode ::DetachSubtree() {
  if (this->pParent == NULL) {
    // nothing needs to be done, since we are trying to sepearet the WHOLE tree
    // it does not make sense....
    return;
  }

  // detach this node (and its descendents) from the rest of the tree
  // note this include free up the current parent
  // this function needs to mantain the coherance of the other tree

  // First seperate the current parent
  RBTNode *pOther = this->pParent->GetLeftChild();
  if (this->IsLeftChild() == true) {
    pOther = this->pParent->GetRightChild();
  }
  pOther->SetParent(this->pParent->GetParent());
  if (this->pParent->GetParent() != NULL) {
    if (this->pParent->IsLeftChild() == true) {
      this->pParent->GetParent()->SetLeftChild(pOther);
    } else {
      this->pParent->GetParent()->SetRightChild(pOther);
    }
  }
  this->pParent->SetLeftChild(NULL);
  this->pParent->SetRightChild(NULL);
  delete this->pParent;

  // need to readjust the tree since the remainig tree may have problem
  // with left/right ordering
  pOther->AdjustLRChildUpwards();

  // Finally set the current node's par to emtpy (meaning detached)
  this->pParent = NULL;
}

RBTNode *RBTNode ::AttachSubtree(RBTNode *pSib) {
  YW_ASSERT_INFO(pSib != NULL, "Fail 2.0");

  // reattach the subtree with its sibling
  // we need to create a new node (which will be returned)
  // this new node could be the new root
  bool fLeftOfSib = true;
  if (this->GetMinLeaveId() > pSib->GetMinLeaveId()) {
    fLeftOfSib = false;
  }
  // cout << "psib = " << (int) pSib << endl;
  // save the original par of psib
  RBTNode *pParSib = pSib->GetParent();

  RBTNode *pPar;
  if (fLeftOfSib == true) {
    pPar = new RBTNode(this, pSib);
  } else {
    pPar = new RBTNode(pSib, this);
  }
  pPar->SetParent(pParSib);
  // cout << "After set parent, create a new node = " << (int) pPar << endl;
  if (pParSib != NULL) {
    if (pParSib->GetLeftChild() == pSib) {
      // cout << "set " << (int) pParSib << " left child to " << (int) pPar <<
      // endl;
      pParSib->SetLeftChild(pPar);
    } else {
      // cout <<"set " << (int) pParSib << " right child to " << (int) pPar <<
      // endl;
      pParSib->SetRightChild(pPar);
    }
  }

  // make sure tree is in right topology
  AdjustLRChildUpwards();

  // cout << "exit from attachsubtree..\n";
  return pPar;
}

RBTNode *RBTNode ::FindLeaf(int lvidParam, int &ponid) {
  // IMPORTANT, in traversal,
  // assume post-order search, and return the how many nodes visited so far
  // Note, ponid should be initialized upon entry (to -1)

  if (IsLeaf() == false) {
    RBTNode *plv = pLeft->FindLeaf(lvidParam, ponid);
    if (plv != NULL) {
      return plv;
    }
    plv = pRight->FindLeaf(lvidParam, ponid);
    if (plv != NULL) {
      return plv;
    }
  }
  // otherwise, increment counter
  ponid++;
  if (IsLeaf() == true) {
    // cout << "visiting leaf = " << this->lvid << ", to search for " <<
    // lvidParam << endl;
    if (this->lvid == lvidParam) {
      return this;
    } else {
      return NULL;
    }
  }
  return NULL;
}

bool RBTNode ::RemoveLeafSelf() {
  // only remove self if it is a leaf
  if (IsLeaf() == false) {
    return false;
  }
  // remove this node
  if (this->pParent != NULL) {
    // need to rearrange the tree to ensure binary shape
    RBTNode *pOther = this->pParent->GetLeftChild();
    if (IsLeftChild() == true) {
      // cout << "Switch to the right\n";
      pOther = this->pParent->GetRightChild();
    }
    // skip the parent
    pOther->SetParent(this->pParent->GetParent());
    // cout << "after getparent\n";
    if (this->pParent->GetParent() != NULL) {
      // cout << "Still need to set parent's parent\n";
      // also need to ensure the proper pointer
      if (pParent->IsLeftChild() == true) {
        pParent->GetParent()->SetLeftChild(pOther);
      } else {
        pParent->GetParent()->SetRightChild(pOther);
      }
    }
    // cout << "delete the old parent\n";
    // free up the parent
    pParent->SetLeftChild(NULL);
    pParent->SetRightChild(NULL);
    delete this->pParent;
    // delete this;

    // make sure the left is ALWAYS smaller than RIGHT
    // BUT SINCE WE ARE REMOVING IN DESCENDING ORDER
    // so it does not matter here. But need to be fixed
    // TBD

  } else {
    // delete this;
  }

  // cout << "done\n";
  return true;
}

// access
int RBTNode ::GetMinLeaveId() {
  YW_ASSERT_INFO(IsLeaf() == true || (pLeft != NULL && pRight != NULL),
                 "Children wrong.");
  if (IsLeaf() == true) {
    return GetLeafId();
  } else {
    int lid = pLeft->GetMinLeaveId();
    int rid = pRight->GetMinLeaveId();
    if (lid < rid) {
      return lid;
    } else {
      return rid;
    }
  }
}

RBTNode *RBTNode ::GetLeftMostChild() {
  RBTNode *pcur = this;
  while (pcur->IsLeaf() == false) {
    pcur = pcur->GetLeftChild();
  }
  return pcur;
}

RBTNode *RBTNode ::GetSibling() {
  if (GetParent() == NULL) {
    return NULL;
  } else {
    if (IsLeftChild() == true) {
      return GetParent()->GetRightChild();
    } else {
      return GetParent()->GetLeftChild();
    }
  }
}

bool RBTNode ::IsLeaf() const { return pLeft == NULL && pRight == NULL; }

int RBTNode ::GetNumLeavesUnder() {
  // cout << "current node = " << (int) this << endl;
  YW_ASSERT_INFO(IsLeaf() == true || (pLeft != NULL && pRight != NULL),
                 "Children wrong.");
  if (IsLeaf() == true) {
    return 1;
  } else {
    return pLeft->GetNumLeavesUnder() + pRight->GetNumLeavesUnder();
  }
}

void RBTNode ::GetLeaves(set<int> &lvs) {
  // cout << "Get leaves so far for node = " << (int) this << ":  ";
  // DumpIntSet( lvs );
  YW_ASSERT_INFO(IsLeaf() == true || (pLeft != NULL && pRight != NULL),
                 "Children wrong.");
  if (IsLeaf() == true) {
    lvs.insert(this->lvid);
  } else {
    pLeft->GetLeaves(lvs);
    pRight->GetLeaves(lvs);
  }
}

bool RBTNode ::IsLeftChild() {
  // if it has no parent, consider left
  if (this->pParent == NULL) {
    return true;
  }
  if (this->pParent->GetLeftChild() == this) {
    return true;
  } else {
    return false;
  }
}

// memory. free recursively
void RBTNode ::Clear() {
  // NOTE: the current node is not deleted!!!!
  // recursively delete
  if (pLeft != NULL) {
    pLeft->Clear();
    delete pLeft;
    pLeft = NULL;
  }
  if (pRight != NULL) {
    pRight->Clear();
    delete pRight;
    pRight = NULL;
  }

  // delete this;
}

void RBTNode ::AdjustLRChildUpwards() {
  // this function re-adjust the left/right subtrees, starting
  // from the current node, and upwards the tree
  // This is because when something is removed, we have to
  // make sure the tree topology is still what is like before:
  // the left subtree must have its min-leaf lower than right
  // subtree
  RBTNode *pcur = this;
  while (pcur != NULL) {
    //
    if (pcur->IsLeaf() == false && pcur->GetLeftChild()->GetMinLeaveId() >
                                       pcur->GetRightChild()->GetMinLeaveId()) {
      // switch it
      RBTNode *ptmp = pcur->GetLeftChild();
      pcur->SetLeftChild(pcur->GetRightChild());
      pcur->SetRightChild(ptmp);
    }

    // trace upwards
    pcur = pcur->GetParent();
  }
}

void RBTNode ::Dump() const {
  // print leaf only
  // this is simply do a post-order traversal
  if (IsLeaf() == true) {
    cout << " " << this->lvid;
    if (GetHeight() >= 0) {
      cout << "[" << GetHeight() << "]";
    }
    cout << " ";
  } else {
    cout << "( ";
    this->GetLeftChild()->Dump();
    this->GetRightChild()->Dump();
    cout << " )";
    if (GetHeight() >= 0) {
      cout << "[" << GetHeight() << "]";
    }
    cout << " ";
  }
}

string RBTNode ::GetNewick() const {
  // if leaf, fill in the leaf id
  if (IsLeaf() == true) {
    char buf[100];
    sprintf(buf, "%d", this->lvid);
    return string(buf);
  } else {
    string s1 = this->GetLeftChild()->GetNewick();
    string s2 = this->GetRightChild()->GetNewick();
    return string("(") + s1 + string(",") + s2 + string(")");
  }
}

void RBTNode ::AddSiblingToLeaf(int lvid) {
  // add a sibling to the current node, which must be a leaf
  YW_ASSERT_INFO(IsLeaf() == true, "Can not add to a non-leaf node");

  // create a new node
  RBTNode *pnode = new RBTNode(lvid);

  // add it
  pnode->AttachSubtree(this);

  // create a new node
  // but remeber the parent first
  //	RBTNode *ppar = GetParent();
  //	bool fLeftChild = IsLeftChild();
  //	YW_ASSERT_INFO( ppar != NULL, "Can not be NULL" );
  //	RBTNode *pinternal = AddSibling( lvid );
  // setup connection
  //	pinternal->SetParent( ppar );
  //	if( fLeftChild == true )
  //	{
  //		// add to the left
  //		ppar->SetLeftChild( pinternal );
  //	}
  //	else
  //	{
  //		ppar->SetRightChild( pinternal );
  //	}
}

void RBTNode ::OutputNodeGML(ofstream &outFile) {
  outFile << "node [\n";
  char name[100];
  // the name is equal to it
  if (IsLeaf() == true) {
    name[0] = 'v';
    sprintf(&name[1], "%d", GetLeafId());
  } else {
    name[0] = ' ';
    name[1] = '\0';
  }
  outFile << "id " << GetLeafId() << endl;
  outFile << "label ";
  OutputQuotedString(outFile, name);
  outFile << endl;
  outFile << "defaultAtrribute   1\n";
  outFile << "]\n";
  // cout << "Output one node: id = " << GetId() << "\n";
  // handle the children
  if (IsLeaf() == false) {
    GetLeftChild()->OutputNodeGML(outFile);
    GetRightChild()->OutputNodeGML(outFile);
  }
}

void RBTNode ::OutputEdgeGML(ofstream &outFile) {
  char name[100];
  int id1 = GetLeafId();
  if (IsLeaf() == false) {
    for (int i = 0; i < 2; ++i) {
      int id2 = GetLeftChild()->GetLeafId();
      if (i == 1) {
        id2 = GetRightChild()->GetLeafId();
      }

      name[0] = ' ';
      name[1] = '\0';
      //		sprintf(&name[1], "%d-%d", id1, id2 );
      // cout << "Output one edge: " << id1 << ", " << id2 << endl;

      outFile << "edge [\n";
      outFile << "source " << id1 << endl;
      outFile << "target  " << id2 << endl;
      outFile << "label ";
      // cout << "edge label = " << name << endl;
      OutputQuotedString(outFile, name);
      outFile << "\n";
      outFile << "]\n";
    }
  }
  // handle the children
  if (IsLeaf() == false) {
    GetLeftChild()->OutputEdgeGML(outFile);
    GetRightChild()->OutputEdgeGML(outFile);
  }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// different ways of initializing a tree
// it can be by a supplied id
RBT ::RBT(int numLeaves, RBT_ID tid) {
  Init();

  // save the id
  this->numLeaves = numLeaves;
  this->tid = tid;
  YW_ASSERT_INFO(numLeaves >= 3, "Too few leaves");

  // construct by the tid
  ReconstructById(tid);
}

RBT ::RBT(const RBT &rhs) {
  this->numLeaves = rhs.numLeaves;
  this->tid = rhs.tid;
  this->pRoot = rhs.pRoot->CopySubTree();
}

RBT &RBT ::operator=(const RBT &rhs) {
  // get rid of current
  if (this->pRoot != NULL) {
    delete this->pRoot;
    this->pRoot = NULL;
  }

  this->numLeaves = rhs.numLeaves;
  this->tid = rhs.tid;
  this->pRoot = rhs.pRoot->CopySubTree();

  return *this;
}

RBT ::RBT(int numLeaves, const vector<int> &listNodeLabels,
          const vector<int> &listParentNodePos,
          const vector<double> &listEdgeDist) {
  this->numLeaves = numLeaves;
  this->tid = -1; // in this mode, we do not care about tid
  this->pRoot = NULL;
  // construct by the tid
  ReconstructByPlainDesc(listNodeLabels, listParentNodePos, listEdgeDist);
}

RBT ::~RBT() {
  // cout << "INside destructor\n";
  // cout << "number of leaves = " << pRoot->GetNumLeavesUnder() << endl;
  this->pRoot->Clear();
  delete pRoot;
  pRoot = NULL;
  // cout << "done with one destructor\n";
}

// ID functions
RBT_ID RBT ::GetId() {
  if (tid >= 0) {
    // return the cached one
    return tid;
  }
  // get it
  this->tid = MapToId(); // indicate it is invalid
  return this->tid;
}

void RBT ::OutputGML(const char *fileName) {
  // Now output a file in GML format
  // First create a new name
  string name = fileName;
  // cout << "num edges = " << listEdges.size() << endl;

  DEBUG("FileName=");
  DEBUG(name);
  DEBUG("\n");
  // Now open file to write out
  ofstream outFile(name.c_str());

  // First output some header info
  outFile << "graph [\n";
  outFile << "comment ";
  OutputQuotedString(outFile, "Automatically generated by Graphing tool");
  outFile << "\ndirected  1\n";
  outFile << "id  1\n";
  outFile << "label ";
  OutputQuotedString(outFile, "To be more meaningful later....\n");
  // cout << "Here we go\n";
  // Now output all the vertices by simply calling through root node
  pRoot->OutputNodeGML(outFile);

  // Now output all the edges by calling through the root
  pRoot->OutputEdgeGML(outFile);

  // Finally quite after closing file
  outFile << "\n]\n";
  outFile.close();
}

// splits functions
bool RBT ::IsSplitContained(const set<int> &split) {
  // simply check the map
  if (mapSplitsInTree.size() == 0) {
    // Need to figure out splits
    RetrieveSplits();
  }
  return mapSplitsInTree.find(split) != mapSplitsInTree.end();
}

void RBT ::GetAllSplits(vector<set<int> > &listSplits) {
  if (mapSplitsInTree.size() == 0) {
    // Need to figure out splits
    RetrieveSplits();
  }

  listSplits.clear();
  for (map<set<int>, bool>::iterator it = mapSplitsInTree.begin();
       it != mapSplitsInTree.end(); ++it) {
    // put it
    listSplits.push_back(it->first);
  }
}

// SPR function
void RBT ::FindSPRDistOneNgbrs(set<int> &ngbrIds) {
  // Double loop: first try every subtree of the original
  // then try to attach it to each of the original node
  // note, we do not want to re-generate trees many times
  // so we need to re-attach the detached subtrees each time we need
  RBT treeOpt(*this);

  TraversRecord tr;
  treeOpt.InitPostorderTranvers(tr);
  while (true) {
    RBTNode *pCurNode = tr.pCurNode;
    // cout << "Outer loop pcurnode = " << (int)pCurNode << ", lvid = " <<
    // pCurNode->GetLeafId() << endl;
    if (pCurNode->GetParent() == NULL) {
      // do not do the whole tree to remove, that is not valid
      break;
    }

    // remember the sibling so we can re-attach it at the end
    RBTNode *pSib = pCurNode->GetParent()->GetLeftChild();
    if (pSib == pCurNode) {
      pSib = pCurNode->GetParent()->GetRightChild();
    }

    // now detach the subtree
    // need to handle the special case when the root is removed
    if (pCurNode->GetParent()->GetParent() == NULL) {
      treeOpt.pRoot = pSib;
    }
    pCurNode->DetachSubtree();
    // set<int> clvs;
    // pCurNode->GetLeaves( clvs );
    // cout << "Current subtree has leafs = ";
    // DumpIntSet( clvs );
    // set<int> rlvs;
    // treeOpt.pRoot->GetLeaves( rlvs );
    // cout << "Remaing tree has leafs = ";
    // DumpIntSet( rlvs );
    // cout << "Current subtree = ";
    // treeOpt.Dump();

    // now do another search
    TraversRecord tr2;
    treeOpt.InitPostorderTranvers(tr2);
    while (true) {
      // set<int> rlvs3;
      // treeOpt.pRoot->GetLeaves( rlvs3 );
      // cout << "During inner loop start, tree has leafs = ";
      // DumpIntSet( rlvs3 );
      // cout << "During internal loop, subtree = ";
      // treeOpt.Dump();

      // cout << "Consider inner node = " << (int)tr2.pCurNode << ", leaf id = "
      // << tr2.pCurNode->GetLeafId() << endl;
      // try to re-attach to the node
      RBTNode *pNewPar = pCurNode->AttachSubtree(tr2.pCurNode);
      if (tr2.pCurNode == treeOpt.pRoot) {
        // we created a new root
        treeOpt.pRoot = pNewPar;
      }

      // get a maped id
      ngbrIds.insert(treeOpt.MapToId());
      // cout << "The SPR transformed subtree = ";
      // treeOpt.Dump();

      // now we need to detach the node again
      if (pCurNode->GetParent()->IsRoot() == true) {
        // when root is removed, we have to re-adjust the root
        treeOpt.pRoot = tr2.pCurNode;
      }
      pCurNode->DetachSubtree();

      // move to next
      if (treeOpt.NextPostorderTranvers(tr2) == false) {
        break;
      }
    }
    // cout << "Now attach the current subtree...\n";
    // now re-attach the node
    RBTNode *pnode = pCurNode->AttachSubtree(pSib);
    if (treeOpt.pRoot == pSib) {
      // cout << "readjust root ...\n";
      // we need to update the root again
      treeOpt.pRoot = pnode;
    }
    // set<int> rlvs2;
    // treeOpt.pRoot->GetLeaves( rlvs2 );
    // cout << "After reattaching at the end of one round, tree has leafs = ";
    // DumpIntSet( rlvs2 );
    // cout << "After re-attaching the subtree = ";
    // treeOpt.Dump();

    // move to next
    if (treeOpt.NextPostorderTranvers(tr) == false) {
      break;
    }
  }

#if 0
	set<RBT> ngbrTrees;
	FindSPRDistOneNgbrs(ngbrTrees);
	for( set<RBT> :: iterator it = ngbrTrees.begin(); it != ngbrTrees.end(); ++it )
	{
		RBT tr = *it;
		ngbrIds.insert( tr.MapToId()  );
	}
#endif
  // get rid of the same tree
  ngbrIds.erase(GetId());
}

void RBT ::FindSPRDistOneNgbrs(vector<RBT *> &ngbrTrees) {
  // Double loop: first try every subtree of the original
  // then try to attach it to each of the original node
  // note, we do not want to re-generate trees many times
  // so we need to re-attach the detached subtrees each time we need
  RBT treeOpt(*this);
  // cout << "RBT: find SPR ngbr: current tree: " << treeOpt.GetNewick() <<
  // endl;

  TraversRecord tr;
  treeOpt.InitPostorderTranvers(tr);
  while (true) {
    RBTNode *pCurNode = tr.pCurNode;
    // cout << "Outer loop pcurnode = " << (int)pCurNode << ", lvid = " <<
    // pCurNode->GetLeafId() << endl;
    if (pCurNode->GetParent() == NULL) {
      // do not do the whole tree to remove, that is not valid
      break;
    }

    // remember the sibling so we can re-attach it at the end
    RBTNode *pSib = pCurNode->GetParent()->GetLeftChild();
    if (pSib == pCurNode) {
      pSib = pCurNode->GetParent()->GetRightChild();
    }

    // now detach the subtree
    // need to handle the special case when the root is removed
    if (pCurNode->GetParent()->GetParent() == NULL) {
      treeOpt.pRoot = pSib;
    }
    pCurNode->DetachSubtree();
    // set<int> clvs;
    // pCurNode->GetLeaves( clvs );
    // cout << "Current subtree has leafs = ";
    // DumpIntSet( clvs );
    // set<int> rlvs;
    // treeOpt.pRoot->GetLeaves( rlvs );
    // cout << "Remaing tree has leafs = ";
    // DumpIntSet( rlvs );
    // cout << "Current subtree = ";
    // treeOpt.Dump();

    // now do another search
    TraversRecord tr2;
    treeOpt.InitPostorderTranvers(tr2);
    while (true) {
      // set<int> rlvs3;
      // treeOpt.pRoot->GetLeaves( rlvs3 );
      // cout << "During inner loop start, tree has leafs = ";
      // DumpIntSet( rlvs3 );
      // cout << "During internal loop, subtree = ";
      // treeOpt.Dump();

      // cout << "Consider inner node = " << (int)tr2.pCurNode << ", leaf id = "
      // << tr2.pCurNode->GetLeafId() << endl;
      // try to re-attach to the node
      RBTNode *pNewPar = pCurNode->AttachSubtree(tr2.pCurNode);
      if (tr2.pCurNode == treeOpt.pRoot) {
        // we created a new root
        treeOpt.pRoot = pNewPar;
      }

      // get a maped id
      // Create a new tree and store it
      RBT *pRbtStore = new RBT(treeOpt);
      ngbrTrees.push_back(pRbtStore);
      // ngbrIds.insert( treeOpt.MapToId()  );
      // cout << "The SPR transformed subtree = ";
      // cout << pRbtStore->GetNewick() << endl;
      // treeOpt.Dump();

      // now we need to detach the node again
      if (pCurNode->GetParent()->IsRoot() == true) {
        // when root is removed, we have to re-adjust the root
        treeOpt.pRoot = tr2.pCurNode;
      }
      pCurNode->DetachSubtree();

      // move to next
      if (treeOpt.NextPostorderTranvers(tr2) == false) {
        break;
      }
    }
    // cout << "Now attach the current subtree...\n";
    // now re-attach the node
    RBTNode *pnode = pCurNode->AttachSubtree(pSib);
    if (treeOpt.pRoot == pSib) {
      // cout << "readjust root ...\n";
      // we need to update the root again
      treeOpt.pRoot = pnode;
    }
    // set<int> rlvs2;
    // treeOpt.pRoot->GetLeaves( rlvs2 );
    // cout << "After reattaching at the end of one round, tree has leafs = ";
    // DumpIntSet( rlvs2 );
    // cout << "After re-attaching the subtree = ";
    // treeOpt.Dump();

    // move to next
    if (treeOpt.NextPostorderTranvers(tr) == false) {
      break;
    }
  }
}

void RBT ::FindSPRDistOneNgbrsRestricted(vector<RBT *> &ngbrTrees,
                                         const vector<RBT *> &ConstraintTrees) {
  // this is slightly different from previous tree in that
  // we want to narrow down on the number of ngbrs to test, thus
  // we want to find more promising ngbrs. In particular,
  // we want to ensure the source branch has a split
  // that is at least one of the constraint trees
  // because the source branch will continue to be one of the splits after
  // transform also, the destination, after merging, the destination new split
  // need to be in one of the constraint tree
  RBT treeOpt(*this);
  int nExcluded = 0;

  TraversRecord tr;
  treeOpt.InitPostorderTranvers(tr);
  while (true) {
    RBTNode *pCurNode = tr.pCurNode;
    // cout << "Outer loop pcurnode = " << (int)pCurNode << ", lvid = " <<
    // pCurNode->GetLeafId() << endl;
    if (pCurNode->GetParent() == NULL) {
      // do not do the whole tree to remove, that is not valid
      break;
    }

    // make sure its leaves are under one of the constriant tree split
    set<int> lvids;
    pCurNode->GetLeaves(lvids);
    // make complmenet if we need
    if (lvids.find(0) == lvids.end()) {
      set<int> tmpset;
      PopulateSetWithInterval(tmpset, 0, this->numLeaves - 1);
      SubtractSets(tmpset, lvids);
      lvids = tmpset;
    }
    bool fContainsrc = false;
    for (int ii = 0; ii < (int)ConstraintTrees.size(); ++ii) {
      RBT *pt = ConstraintTrees[ii];
      YW_ASSERT_INFO(pt != NULL, "wrong");
      if (pt->IsSplitContained(lvids) == true) {
        fContainsrc = true;
        break;
      }
    }
    if (fContainsrc == false) {
      nExcluded++;
    }

    if (fContainsrc == true) {

      // remember the sibling so we can re-attach it at the end
      RBTNode *pSib = pCurNode->GetParent()->GetLeftChild();
      if (pSib == pCurNode) {
        pSib = pCurNode->GetParent()->GetRightChild();
      }

      // now detach the subtree
      // need to handle the special case when the root is removed
      if (pCurNode->GetParent()->GetParent() == NULL) {
        treeOpt.pRoot = pSib;
      }
      pCurNode->DetachSubtree();
      // set<int> clvs;
      // pCurNode->GetLeaves( clvs );
      // cout << "Current subtree has leafs = ";
      // DumpIntSet( clvs );
      // set<int> rlvs;
      // treeOpt.pRoot->GetLeaves( rlvs );
      // cout << "Remaing tree has leafs = ";
      // DumpIntSet( rlvs );
      // cout << "Current subtree = ";
      // treeOpt.Dump();

      // now do another search
      TraversRecord tr2;
      treeOpt.InitPostorderTranvers(tr2);
      while (true) {
        // set<int> rlvs3;
        // treeOpt.pRoot->GetLeaves( rlvs3 );
        // cout << "During inner loop start, tree has leafs = ";
        // DumpIntSet( rlvs3 );
        // cout << "During internal loop, subtree = ";
        // treeOpt.Dump();

        // cout << "Consider inner node = " << (int)tr2.pCurNode << ", leaf id =
        // " << tr2.pCurNode->GetLeafId() << endl;
        // try to re-attach to the node
        RBTNode *pNewPar = pCurNode->AttachSubtree(tr2.pCurNode);
        if (tr2.pCurNode == treeOpt.pRoot) {
          // we created a new root
          treeOpt.pRoot = pNewPar;
        }

        // is pNewPar has a split that exists in one of constraint tree?
        set<int> lvids2;
        pNewPar->GetLeaves(lvids2);
        // make complmenet if we need
        if (lvids2.find(0) == lvids2.end()) {
          set<int> tmpset;
          PopulateSetWithInterval(tmpset, 0, this->numLeaves - 1);
          SubtractSets(tmpset, lvids2);
          lvids2 = tmpset;
        }
        bool fContainsrc2 = false;
        for (int ii = 0; ii < (int)ConstraintTrees.size(); ++ii) {
          RBT *pt = ConstraintTrees[ii];
          YW_ASSERT_INFO(pt != NULL, "wrong");
          if (pt->IsSplitContained(lvids2) == true) {
            fContainsrc2 = true;
            break;
          }
        }
        if (fContainsrc2 == true) {
          // get a maped id
          // Create a new tree and store it
          RBT *pRbtStore = new RBT(treeOpt);
          ngbrTrees.push_back(pRbtStore);
          // ngbrIds.insert( treeOpt.MapToId()  );
          // cout << "The SPR transformed subtree = ";
          // treeOpt.Dump();
        }

        // now we need to detach the node again
        if (pCurNode->GetParent()->IsRoot() == true) {
          // when root is removed, we have to re-adjust the root
          treeOpt.pRoot = tr2.pCurNode;
        }
        pCurNode->DetachSubtree();

        // move to next
        if (treeOpt.NextPostorderTranvers(tr2) == false) {
          break;
        }
      }
      // cout << "Now attach the current subtree...\n";
      // now re-attach the node
      RBTNode *pnode = pCurNode->AttachSubtree(pSib);
      if (treeOpt.pRoot == pSib) {
        // cout << "readjust root ...\n";
        // we need to update the root again
        treeOpt.pRoot = pnode;
      }
      // set<int> rlvs2;
      // treeOpt.pRoot->GetLeaves( rlvs2 );
      // cout << "After reattaching at the end of one round, tree has leafs = ";
      // DumpIntSet( rlvs2 );
      // cout << "After re-attaching the subtree = ";
      // treeOpt.Dump();
    }

    // move to next
    if (treeOpt.NextPostorderTranvers(tr) == false) {
      break;
    }
  }

  cout << "excluded num = " << nExcluded << endl;
}

// is tree SPR away from this
bool RBT ::IsOneSPRAway(const RBT &rbt) const {
  // testing whether it is one SPR away
  // Simply try to morph the current tree t
  // Double loop: first try every subtree of the original
  // then try to attach it to each of the original node
  // note, we do not want to re-generate trees many times
  // so we need to re-attach the detached subtrees each time we need
  // BUT, to make process fast, we need to reduce the tree as much as we can
  //
  RBT treeOpt(*this);
  RBT treeCmp(rbt);

  // reduce the two trees
  Consolidate(treeOpt, treeCmp);
  // cout <<"After consolidation, trees are: \n";
  // treeOpt.Dump();
  // treeCmp.Dump();

  // first make an list of maps to nodes at tips
  treeOpt.CollectTips();
  treeCmp.CollectTips();
  vector<RBTNode *> listTips1;
  treeOpt.GetAllTips(listTips1);
  // cout << "Find tip num = " << listTips1.size() << endl;
  // store all pair of nodes s.t. it only appears in treeOpt
  // In fact, if the preprocessing step is correct,
  // a cherry (a pair of nodes) appears in treeA can NOT appear in treeB
  map<pair<RBTNode *, RBTNode *>, bool> mapCherry1;
  for (int i = 0; i < (int)listTips1.size(); ++i) {
    // cout << "Processing tip = " << listTips1[i]->GetLeafId() << endl;
    // get its sibling
    RBTNode *pSib = listTips1[i]->GetSibling();
    if (pSib->IsLeaf() == true) {
      // cout << "Sibling is a LEAF...\n";
      pair<RBTNode *, RBTNode *> pp;
      // get rid of order
      if ((long)pSib > (long)listTips1[i]) {
        pp.first = listTips1[i];
        pp.second = pSib;
      } else {
        pp.second = listTips1[i];
        pp.first = pSib;
      }
      mapCherry1.insert(
          map<pair<RBTNode *, RBTNode *>, bool>::value_type(pp, true));

      // make sure preprocessing is correct
      // by checking the situation at the other tree
      // the same pair can NOT appear
      RBTNode *pOther1 = treeCmp.GetTip(pp.first->GetLeafId());
      RBTNode *pOtherSib = pOther1->GetSibling();
      RBTNode *pOther2 = treeCmp.GetTip(pp.second->GetLeafId());
      YW_ASSERT_INFO(pOtherSib != pOther2, "Tree preprocessing wrong");
    }
  }
  // if there is more than 2 pair left, we are done
  if (mapCherry1.size() >= 3) {
    //
    return false;
  }
  YW_ASSERT_INFO(mapCherry1.size() > 0 && mapCherry1.size() < 3,
                 "Wrong: cherry number can not be empty");
  //  In this case, pick one pair (say the first), and perform one SPR to get a
  //  proper subset
  // collect the list of leaf edges to try
  // vector< RBTNode *> listLeafToBePruned, listRegraftDest;
  // for( map< pair<RBTNode *, RBTNode *>, bool > :: iterator it =
  // mapCherry1.begin(); it != mapCherry1.end(); ++it )
  //{
  //	listLeafToBePruned.push_back( it->first.first );
  //	listLeafToBePruned.push_back( it->first.second );
  //}
  // also figure out the destination it has to be
  // for(int i=0; i<(int)listLeafToBePruned.size();++i)
  //{
  //
  //}

  // first, if there is only one pair of tips, then the tree must be like a comb
  RBTNode *pLeaf1 = NULL;
  RBTNode *pLeaf2 = NULL;
  RBTNode *pLeaf3 = NULL;
  RBTNode *pLeaf4 = NULL;
  map<pair<RBTNode *, RBTNode *>, bool>::iterator it = mapCherry1.begin();
  pLeaf1 = it->first.first;
  pLeaf2 = it->first.second;
  it++;
  if (it != mapCherry1.end()) {
    pLeaf3 = it->first.first;
    pLeaf4 = it->first.second;
  }

  // now start real comparasion
  TraversRecord tr;
  treeOpt.InitPostorderTranvers(tr);
  while (true) {
    RBTNode *pCurNode = tr.pCurNode;
    // cout << "Outer loop pcurnode = " << (int)pCurNode << ", lvid = " <<
    // pCurNode->GetLeafId() << endl;
    if (pCurNode->GetParent() == NULL) {
      // do not do the whole tree to remove, that is not valid
      break;
    }

    // remember the sibling so we can re-attach it at the end
    RBTNode *pSib = pCurNode->GetParent()->GetLeftChild();
    if (pSib == pCurNode) {
      pSib = pCurNode->GetParent()->GetRightChild();
    }

    // now detach the subtree
    // need to handle the special case when the root is removed
    if (pCurNode->GetParent()->GetParent() == NULL) {
      treeOpt.pRoot = pSib;
    }
    pCurNode->DetachSubtree();
    // set<int> clvs;
    // pCurNode->GetLeaves( clvs );
    // cout << "Current subtree has leafs = ";
    // DumpIntSet( clvs );
    // set<int> rlvs;
    // treeOpt.pRoot->GetLeaves( rlvs );
    // cout << "Remaing tree has leafs = ";
    // DumpIntSet( rlvs );
    // cout << "Current subtree = ";
    // treeOpt.Dump();

    // now do another search
    TraversRecord tr2;
    treeOpt.InitPostorderTranvers(tr2);
    while (true) {
      // set<int> rlvs3;
      // treeOpt.pRoot->GetLeaves( rlvs3 );
      // cout << "During inner loop start, tree has leafs = ";
      // DumpIntSet( rlvs3 );
      // cout << "During internal loop, subtree = ";
      // treeOpt.Dump();

      // make sure this node is what we need:
      // (1) must be a leaf
      if (tr2.pCurNode->IsLeaf() == true &&
          ((mapCherry1.size() == 1 &&
            (pCurNode == pLeaf1 || pCurNode == pLeaf2 ||
             tr2.pCurNode == pLeaf1 || tr2.pCurNode == pLeaf2)) ||
           (mapCherry1.size() == 2 &&
            (((pCurNode == pLeaf1 || pCurNode == pLeaf2) &&
              (tr2.pCurNode == pLeaf3 || tr2.pCurNode == pLeaf4)) ||
             ((pCurNode == pLeaf3 || pCurNode == pLeaf4) &&
              (tr2.pCurNode == pLeaf1 || tr2.pCurNode == pLeaf2)))))) {

        // cout << "Consider inner node = " << (int)tr2.pCurNode << ", leaf id =
        // " << tr2.pCurNode->GetLeafId() << endl;
        // try to re-attach to the node
        RBTNode *pNewPar = pCurNode->AttachSubtree(tr2.pCurNode);
        if (tr2.pCurNode == treeOpt.pRoot) {
          // we created a new root
          treeOpt.pRoot = pNewPar;
        }

        // Test whether the morphed tree is the SAME as the other
        if (treeOpt.IsSame(treeCmp) == true) {
          // find it
          return true;
          // cout << "The SPR transformed subtree = ";
          // treeOpt.Dump();
        }

        // now we need to detach the node again
        if (pCurNode->GetParent()->IsRoot() == true) {
          // when root is removed, we have to re-adjust the root
          treeOpt.pRoot = tr2.pCurNode;
        }
        pCurNode->DetachSubtree();
      }
      // move to next
      if (treeOpt.NextPostorderTranvers(tr2) == false) {
        break;
      }
    }
    // cout << "Now attach the current subtree...\n";
    // now re-attach the node
    RBTNode *pnode = pCurNode->AttachSubtree(pSib);
    if (treeOpt.pRoot == pSib) {
      // cout << "readjust root ...\n";
      // we need to update the root again
      treeOpt.pRoot = pnode;
    }
    // set<int> rlvs2;
    // treeOpt.pRoot->GetLeaves( rlvs2 );
    // cout << "After reattaching at the end of one round, tree has leafs = ";
    // DumpIntSet( rlvs2 );
    // cout << "After re-attaching the subtree = ";
    // treeOpt.Dump();

    // move to next
    if (treeOpt.NextPostorderTranvers(tr) == false) {
      break;
    }
  }

  // did not find
  return false;

#if 0
	// now start real comparasion
    TraversRecord tr;
    treeOpt.InitPostorderTranvers(tr);
    while( true )
    {
        RBTNode *pCurNode = tr.pCurNode;
//cout << "Outer loop pcurnode = " << (int)pCurNode << ", lvid = " << pCurNode->GetLeafId() << endl;
        if( pCurNode->GetParent() == NULL  )
        {
            // do not do the whole tree to remove, that is not valid
            break;
        }

        // remember the sibling so we can re-attach it at the end
        RBTNode *pSib = pCurNode->GetParent()->GetLeftChild();
        if( pSib == pCurNode )
        {
            pSib = pCurNode->GetParent()->GetRightChild();
        }



        // now detach the subtree
        // need to handle the special case when the root is removed
        if( pCurNode->GetParent()->GetParent() == NULL )
        {
            treeOpt.pRoot = pSib;
        }
        pCurNode->DetachSubtree();
//set<int> clvs;
//pCurNode->GetLeaves( clvs );
//cout << "Current subtree has leafs = ";
//DumpIntSet( clvs );
//set<int> rlvs;
//treeOpt.pRoot->GetLeaves( rlvs );
//cout << "Remaing tree has leafs = ";
//DumpIntSet( rlvs );
//cout << "Current subtree = ";
//treeOpt.Dump();

        // now do another search
        TraversRecord tr2;
        treeOpt.InitPostorderTranvers( tr2 );
        while(true)
        {
//set<int> rlvs3;
//treeOpt.pRoot->GetLeaves( rlvs3 );
//cout << "During inner loop start, tree has leafs = ";
//DumpIntSet( rlvs3 );
//cout << "During internal loop, subtree = ";
//treeOpt.Dump();


//cout << "Consider inner node = " << (int)tr2.pCurNode << ", leaf id = " << tr2.pCurNode->GetLeafId() << endl;
            // try to re-attach to the node
            RBTNode *pNewPar = pCurNode->AttachSubtree(tr2.pCurNode);
            if( tr2.pCurNode == treeOpt.pRoot )
            {
                // we created a new root
                treeOpt.pRoot = pNewPar;
            }

            // Test whether the morphed tree is the SAME as the other
            if (treeOpt.IsSame( treeCmp ) == true  )
			{
				// find it
				return true;
//cout << "The SPR transformed subtree = ";
//treeOpt.Dump();
			}

            // now we need to detach the node again
            if( pCurNode->GetParent()->IsRoot() == true  )
            {
                // when root is removed, we have to re-adjust the root
                treeOpt.pRoot = tr2.pCurNode;
            }
            pCurNode->DetachSubtree();

            // move to next
            if( treeOpt.NextPostorderTranvers(tr2) == false )
            {
                break;
            }

        }
//cout << "Now attach the current subtree...\n";
        // now re-attach the node
        RBTNode *pnode = pCurNode->AttachSubtree( pSib );
        if( treeOpt.pRoot == pSib )
        {
//cout << "readjust root ...\n";
            // we need to update the root again
            treeOpt.pRoot = pnode;
        }
//set<int> rlvs2;
//treeOpt.pRoot->GetLeaves( rlvs2 );
//cout << "After reattaching at the end of one round, tree has leafs = ";
//DumpIntSet( rlvs2 );
//cout << "After re-attaching the subtree = ";
//treeOpt.Dump();

        // move to next
        if( treeOpt.NextPostorderTranvers(tr) == false )
        {
            break;
        }
    }

	// did not find
	return false;
#endif

#if 0
	// testing whether it is one SPR away
	// Simply try to morph the current tree t
    // Double loop: first try every subtree of the original
    // then try to attach it to each of the original node
    // note, we do not want to re-generate trees many times
    // so we need to re-attach the detached subtrees each time we need
    RBT treeOpt(*this);


    TraversRecord tr;
    treeOpt.InitPostorderTranvers(tr);
    while( true )
    {
        RBTNode *pCurNode = tr.pCurNode;
//cout << "Outer loop pcurnode = " << (int)pCurNode << ", lvid = " << pCurNode->GetLeafId() << endl;
        if( pCurNode->GetParent() == NULL  )
        {
            // do not do the whole tree to remove, that is not valid
            break;
        }

        // remember the sibling so we can re-attach it at the end
        RBTNode *pSib = pCurNode->GetParent()->GetLeftChild();
        if( pSib == pCurNode )
        {
            pSib = pCurNode->GetParent()->GetRightChild();
        }



        // now detach the subtree
        // need to handle the special case when the root is removed
        if( pCurNode->GetParent()->GetParent() == NULL )
        {
            treeOpt.pRoot = pSib;
        }
        pCurNode->DetachSubtree();
//set<int> clvs;
//pCurNode->GetLeaves( clvs );
//cout << "Current subtree has leafs = ";
//DumpIntSet( clvs );
//set<int> rlvs;
//treeOpt.pRoot->GetLeaves( rlvs );
//cout << "Remaing tree has leafs = ";
//DumpIntSet( rlvs );
//cout << "Current subtree = ";
//treeOpt.Dump();

        // now do another search
        TraversRecord tr2;
        treeOpt.InitPostorderTranvers( tr2 );
        while(true)
        {
//set<int> rlvs3;
//treeOpt.pRoot->GetLeaves( rlvs3 );
//cout << "During inner loop start, tree has leafs = ";
//DumpIntSet( rlvs3 );
//cout << "During internal loop, subtree = ";
//treeOpt.Dump();


//cout << "Consider inner node = " << (int)tr2.pCurNode << ", leaf id = " << tr2.pCurNode->GetLeafId() << endl;
            // try to re-attach to the node
            RBTNode *pNewPar = pCurNode->AttachSubtree(tr2.pCurNode);
            if( tr2.pCurNode == treeOpt.pRoot )
            {
                // we created a new root
                treeOpt.pRoot = pNewPar;
            }

            // Test whether the morphed tree is the SAME as the other
            if (treeOpt.IsSame( rbt ) == true  )
			{
				// find it
				return true;
//cout << "The SPR transformed subtree = ";
//treeOpt.Dump();
			}

            // now we need to detach the node again
            if( pCurNode->GetParent()->IsRoot() == true  )
            {
                // when root is removed, we have to re-adjust the root
                treeOpt.pRoot = tr2.pCurNode;
            }
            pCurNode->DetachSubtree();

            // move to next
            if( treeOpt.NextPostorderTranvers(tr2) == false )
            {
                break;
            }

        }
//cout << "Now attach the current subtree...\n";
        // now re-attach the node
        RBTNode *pnode = pCurNode->AttachSubtree( pSib );
        if( treeOpt.pRoot == pSib )
        {
//cout << "readjust root ...\n";
            // we need to update the root again
            treeOpt.pRoot = pnode;
        }
//set<int> rlvs2;
//treeOpt.pRoot->GetLeaves( rlvs2 );
//cout << "After reattaching at the end of one round, tree has leafs = ";
//DumpIntSet( rlvs2 );
//cout << "After re-attaching the subtree = ";
//treeOpt.Dump();

        // move to next
        if( treeOpt.NextPostorderTranvers(tr) == false )
        {
            break;
        }
    }

	// did not find
	return false;
#endif
}

// 11/15/07: found an error: sometimes it passed in an invalid tree, then we get
// problems TBD. need to figure out why this is happening This function is to
// reduce two trees, such that the two trees' common parts are removed, only
// different parts are left
void RBT ::Consolidate(RBT &treeOpt, RBT &treeCmp) {
  // cout << "ENTERING consolidate....\n";
  YW_ASSERT_INFO(treeOpt.GetNodesNum() == treeCmp.GetNodesNum(),
                 "Tree must be the same");
  // create a map of leaf nodes ofr cmp tree
  map<int, RBTNode *> mapCmpTreeLeafNodes;
  TraversRecord tr1;
  treeCmp.InitPostorderTranvers(tr1);
  while (true) {
    if (tr1.pCurNode->IsLeaf() == true) {
      mapCmpTreeLeafNodes.insert(map<int, RBTNode *>::value_type(
          tr1.pCurNode->GetLeafId(), tr1.pCurNode));
    }
    if (treeCmp.NextPostorderTranvers(tr1) == false) {
      break;
    }
  }

  // cout << "here1\n";
  // reduce the two trees so that there is shared subtrees in them
  // I do not understand why ONE-PATH does not work. Here just repeat until no
  // nodes can be deleted
  bool fNothingFound = false;
  while (fNothingFound == false) {
    // cout << "Current trees = ";
    // treeOpt.Dump();
    // treeCmp.Dump();
    fNothingFound = true;
    TraversRecord tr;
    treeOpt.InitPostorderTranvers(tr);
    bool fNodeDeleted = false;
    while (true) {
      //
      if (tr.pCurNode->IsLeaf() == true) {
        //
        // if( tr.pCurNode->IsLeftChild() == true )
        // we we start to delete, we only look for left child for now
        YW_ASSERT_INFO(tr.pCurNode->GetParent() != NULL,
                       "Can not be like this");
        RBTNode *psib = tr.pCurNode->GetSibling();
        YW_ASSERT_INFO(psib != NULL, "Wrong1.1.0");
        if (psib->IsLeaf() == true) {
          // now try to get the corresponding leaf in the other tree
          RBTNode *pLeaf1Cmp = mapCmpTreeLeafNodes[tr.pCurNode->GetLeafId()];
          RBTNode *pLeaf2Cmp = mapCmpTreeLeafNodes[psib->GetLeafId()];
          if (pLeaf1Cmp == NULL) {
            // treeOpt.Dump();
            // treeCmp.Dump();
            cout << "This node has been delted: " << tr.pCurNode->GetLeafId()
                 << endl;
          }
          if (pLeaf2Cmp == NULL) {
            // treeOpt.Dump();
            // treeCmp.Dump();
            cout << "This node has been delted: " << psib->GetLeafId() << endl;
          }
          // YW: for now, continue, need to fix it later. 11/15/07
          YW_ASSERT_INFO(pLeaf1Cmp != NULL && pLeaf2Cmp != NULL, "Wrong1.1.1");
          if (pLeaf1Cmp->GetParent() == pLeaf2Cmp->GetParent()) {

            // Good, we find a pair, now we remove the right node
            fNodeDeleted = true;
            fNothingFound = false;
            int sibidCmp = psib->GetLeafId();
            pLeaf2Cmp->RemoveLeafSelf();
            delete pLeaf2Cmp;
            pLeaf2Cmp = NULL;
            mapCmpTreeLeafNodes[sibidCmp] = NULL;

            psib->RemoveLeafSelf();
            delete psib;
            psib = NULL;
            // cout << "Leaf " << sibidCmp << " is deleted\n";
            // if( tr.pCurNode->IsLeftChild() == false )
            // if( tr.pCurNode-> )
            //{
            // Update current to left child
            //
            //	tr.pCurNode = psib;
            //}
          }
        }
      }

      if (fNodeDeleted == true) {
        // give one more chance
        fNodeDeleted = false;
        continue;
      }

      if (treeOpt.NextPostorderTranvers(tr) == false) {
        break;
      }
    }
  }
  // cout << "here2\n";
}

bool RBT ::ReconstructNewick(const string &strNewick) {
  // for now, call internal
  RBTNode *pRootNew = ReconstructNewickInternal(strNewick);
  if (pRootNew == NULL) {
    // fail to build
    return false;
  }
  // update current node
  if (this->pRoot != NULL) {
    this->pRoot->Clear();
    delete pRoot;
  }
  this->pRoot = pRootNew;
  return true;
}

void RBT ::CollectTips() {
  mapTipPtrs.clear();

  //
  TraversRecord tr;
  InitPostorderTranvers(tr);
  while (true) {
    //
    if (tr.pCurNode->IsLeaf() == true) {
      mapTipPtrs.insert(map<int, RBTNode *>::value_type(
          tr.pCurNode->GetLeafId(), tr.pCurNode));
    }

    // continue
    if (NextPostorderTranvers(tr) == false) {
      break;
    }
  }
}

RBTNode *RBT ::GetTip(int id) {
  if (mapTipPtrs.find(id) != mapTipPtrs.end()) {
    return mapTipPtrs[id];
  } else {
    return NULL;
  }
}

void RBT ::GetAllTips(vector<RBTNode *> &tips) {
  for (map<int, RBTNode *>::iterator it = mapTipPtrs.begin();
       it != mapTipPtrs.end(); ++it) {
    tips.push_back(it->second);
  }
}

bool RBT ::AddLeaf(int pos) {
  // make sure this is a good position
  if (pos >= 2 * numLeaves - 1) {
    // bad position
    return false;
  }

  // now add to the leaf
  InternalAddleaf(numLeaves, pos);

  // inc num of leaves
  numLeaves++;

  // clean up
  mapSplitsInTree.clear();
  this->tid = MapToId();
  return true;
}

// compare
int RBT ::Compare(RBT &rhs) {
  // simply find how many splits are common in two trees
  // collect two sets of splits
  vector<set<int> > listSplitsRHS;
  rhs.GetAllSplits(listSplitsRHS);
  set<set<int> > setSplitsRHS;
  for (int i = 0; i < (int)listSplitsRHS.size(); ++i) {
    setSplitsRHS.insert(listSplitsRHS[i]);
  }
  vector<set<int> > listSplits;
  this->GetAllSplits(listSplits);
  int res = 0;
  for (int i = 0; i < (int)listSplits.size(); ++i) {
    if (setSplitsRHS.find(listSplits[i]) != setSplitsRHS.end()) {
      // find oe shared
      res++;
    }
  }
  return res;
}
bool RBT ::IsSameUnrootedTree(RBT &rhs) {
  // simply find how many splits are common in two trees
  // collect two sets of splits
  vector<set<int> > listSplitsRHS;
  rhs.GetAllSplits(listSplitsRHS);
  set<set<int> > setSplitsRHS;
  for (int i = 0; i < (int)listSplitsRHS.size(); ++i) {
    setSplitsRHS.insert(listSplitsRHS[i]);
  }
  vector<set<int> > listSplits;
  this->GetAllSplits(listSplits);
  for (int i = 0; i < (int)listSplits.size(); ++i) {
    if (setSplitsRHS.find(listSplits[i]) == setSplitsRHS.end()) {
      // find oe shared
      return false;
    }
  }
  return true;
}

///////////////////////////////////////////////////////////////////////////////////////

RBTNode *RBT ::ReconstructNewickInternal(const string &strNewick) {
  // Build RBT by a given Newick string
  // NOTE: we assume the tree is in the form of (1,(2,3)) form
  // THAT IS, WE DO NOT ALLOW PRECEEDING SYMBOLS
  // return the constructed root node for the current substring
  // define commonly used symbol in Newick
  // const char cTerm = ';';

  // this function builds recursively subtrees for this part of string
  // First, is this string a leaf or not
  if (strNewick[0] != '(') {
    // Yes, this is a leaf
    int nodeId;
    sscanf(strNewick.c_str(), "%d", &nodeId);
    // cout << "leaf id = " << nodeId << endl;

    // the ID of ms is by convention, one larger (starting from 1)
    // so decrement by one

    RBTNode *pLeaf = new RBTNode(nodeId - 1);
    return pLeaf;
  } else {
    // This is not a leaf
    // so we create underlying level for it
    // TreeNode *pInternal = new TreeNode( invId++  );
    RBTNode *pLeftChild = NULL;
    RBTNode *pRightChild = NULL;
    int lastpos = 1;
    int curpos = 0;
    int parnet = 0; // (: +1, ) -1
    while (true) {
      // cout << "curpos = " << curpos << endl;

      if (curpos >= (int)strNewick.size()) {
        // we are done
        break;
      }

      // keep balance
      if (strNewick[curpos] == '(') {
        parnet++;
      } else if (strNewick[curpos] == ')') {
        parnet--;

        // when parnet = 0, we know we end
        if (parnet == 0) {
          // now adding the last piece
          // create a new node
          int strl = curpos - lastpos;
          string subs = strNewick.substr(lastpos, strl);
          //    cout << "last subs = " << subs << endl;
          pLeftChild = ReconstructNewickInternal(subs);

          // aslo update lastpos
          lastpos = curpos + 1;
        }

      } else if (strNewick[curpos] == ',') {
        // Yes, this is a sepeartor, but we only start to process it when the
        // balance of parenetnis is right
        if (parnet == 1) {
          // create a new node
          int strl = curpos - lastpos;
          string subs = strNewick.substr(lastpos, strl);
          //    cout << "subs = " << subs << endl;
          pRightChild = ReconstructNewickInternal(subs);

          // aslo update lastpos
          lastpos = curpos + 1;
        }
      }

      // now move to next pos
      curpos++;
    }

    YW_ASSERT_INFO(pLeftChild != NULL && pRightChild != NULL, "Children wrong");
    RBTNode *pInternal;
    if (pLeftChild->GetMinLeaveId() < pRightChild->GetMinLeaveId()) {
      pInternal = new RBTNode(pLeftChild, pRightChild);
    } else {
      pInternal = new RBTNode(pRightChild, pLeftChild);
    }
    return pInternal;
  }

  // reconstruct tree by the given Newick format
  // int spos = 0;
  // while( spos < (int) strNewick.size()   )
  //{
  //    if(  strNewick[spos] == cTerm  )
  //    {
  //        break;
  //    }
  //    // Skip things until we find the first (
  //}
}

/////////////////////////////////////////////////////////////////////////////////////
void RBT ::Init() {
  pRoot = NULL;
  tid = -1; // not initialized
  numLeaves = 0;
}

void RBT ::ReconstructById(RBT_ID tid) {
  // cout << "ReconstructById\n";
  // first clear the old tree if any
  if (pRoot != NULL) {
    pRoot->Clear();
    delete pRoot;
    pRoot = NULL;
  }

  vector<int> leavesEdgeIndices(numLeaves);
  leavesEdgeIndices[0] = 0;
  leavesEdgeIndices[1] = 0;

  // reconstruct the tree by its ID
  // first restrive the edge ids
  int idUse = tid;
  for (int lv = numLeaves - 1; lv >= 2; --lv) {
    int base = 2 * lv - 1;
    int eid = idUse % base;
    leavesEdgeIndices[lv] = eid;
    idUse = idUse / base;
  }
  // create a tree with two leaves
  RBTNode *pn0 = new RBTNode(0);
  // cout << "pn0 = " << (int) pn0 << endl;
  RBTNode *pn1 = new RBTNode(1);
  // cout << "pn1 = " << (int) pn1 << endl;
  RBTNode *prn = new RBTNode(pn0, pn1);
  // cout << "prn = " << (int) prn << endl;
  this->pRoot = prn;

  // now start to insert nodes from the third leaf
  for (int lv = 2; lv < numLeaves; ++lv) {
    // cout << "lv = " << lv << ", in construction\n";
    // make sure the index make sense
    int eid = leavesEdgeIndices[lv];
    YW_ASSERT_INFO(eid < 2 * lv - 1, "eid too large");

    InternalAddleaf(lv, eid);

    // cout << "eid = " << eid << endl;
    /*
            // travere the current tree, and stop at the index
            TraversRecord tr;
            InitPostorderTranvers(tr);
            int cureid = 0;
            while(true)
            {
                if( cureid == eid )
                {
                    // find it!
                    break;
                }
                else
                {
                    // continue
                    NextPostorderTranvers(tr);
                }

                if( cureid >= 2*lv-1 )
                {
                    // should not come here
                    YW_ASSERT_INFO(false, "Should not be here");
                    break;
                }
                // update
                cureid ++;
            }
    //cout << "cureid = " << cureid << endl;
            // now add this. Need to consider whether this is the root or not
            if( tr.pCurNode == pRoot )
            {
                RBTNode *pNewRoot = pRoot->AddSibling( lv );
                this->pRoot = pNewRoot;
    //cout << "Update root to " << (int)pNewRoot << endl;
            }
            else
            {
                // then simply add it to parent's proper position
                if( tr.pCurNode->IsLeftChild() == true )
                {
                    tr.pCurNode->GetParent()->AddToLeftEdge(lv);
                }
                else
                {
                    tr.pCurNode->GetParent()->AddToRightEdge(lv);
                }
            }
            */
  }

  // before return, save the clusters
  // TBD
  // YW_ASSERT_INFO(false, "not implemented");
}

// handle insertion of a new leaf
// note: we only allow SEQENTIALLY INSERTION OF LEAVES
bool RBT ::InternalAddleaf(int lvid, int pos) {
  // travere the current tree, and stop at the index
  TraversRecord tr;
  InitPostorderTranvers(tr);
  int cureid = 0;
  while (true) {
    if (cureid == pos) {
      // find it!
      break;
    } else {
      // continue
      NextPostorderTranvers(tr);
    }

    if (cureid >= 2 * lvid - 1) {
      // should not come here
      YW_ASSERT_INFO(false, "Should not be here2");
      break;
    }
    // update
    cureid++;
  }
  // cout << "cureid = " << cureid << endl;
  // now add this. Need to consider whether this is the root or not
  if (tr.pCurNode == pRoot) {
    RBTNode *pNewRoot = pRoot->AddSibling(lvid);
    this->pRoot = pNewRoot;
    // cout << "Update root to " << (int)pNewRoot << endl;
  } else {
    // then simply add it to parent's proper position
    if (tr.pCurNode->IsLeftChild() == true) {
      tr.pCurNode->GetParent()->AddToLeftEdge(lvid);
    } else {
      tr.pCurNode->GetParent()->AddToRightEdge(lvid);
    }
  }
  return true;
}

RBT_ID RBT ::MapToId() {
  // The scheme needs to be carefully worked out
  // We use the enumeration index of the leave as the id base
  // That is, id = [id2, id3, id4, ..., idk]
  // where idi indicates which edge we pick in the RBT when inserting leaf-i
  // We need to choose a way to assign number to (partial-completed)-tree edges
  // we do so by post-order traversal: an edge is assign the POT order to the
  // corresponding node (as the one towards the leaves of the tree)

  YW_ASSERT_INFO(numLeaves >= 3, "Too few leaves");
  // map the tree to an ID
  // we save a vector of indices, which indicates on which edge the split is
  // from
  vector<int> leavesEdgeIndices(numLeaves);
  leavesEdgeIndices[0] = 0;
  leavesEdgeIndices[1] = 0;

  // reconstruct a new tree by copying
  RBT treeNew(*this);
  // cout << "Tree copied. \n";
  // start from third leave
  for (int lv = numLeaves - 1; lv >= 2; --lv) {
    // cout << "lv = " << lv << endl;
    // find out where is this leave
    int ponid = -1;
    RBTNode *pLeaf = treeNew.FindLeaf(lv, ponid);
    YW_ASSERT_INFO(pLeaf != NULL, "Fail in getting a leaf");
    // cout << "ponid = " << ponid << endl;
    if (pLeaf->IsLeftChild() == true) {
      // if LEFt child, then ponid is TRUE
      // so no change here
    } else {
      // cout << "It is right child\n";
      // if it is RIGHT child, then in the original insert,
      // it is put at ponid-1 edge
      ponid--;
    }
    // remmeber this ponid
    leavesEdgeIndices[lv] = ponid;
    // remove this lv
    // here is not very robust, but since we are not deleting the whole thing
    // so it should be OK
    // cout << "leaf id = " << pLeaf->GetLeafId() << endl;

    // update root
    if (pLeaf->GetParent() != NULL && pLeaf->GetParent()->GetParent() == NULL) {
      // cout << "UPdate root\n";
      // in this case, update pRoot
      if (pLeaf->IsLeftChild() == true) {
        treeNew.pRoot = pLeaf->GetParent()->GetRightChild();
      } else {
        // cout << "Get left child\n";
        treeNew.pRoot = pLeaf->GetParent()->GetLeftChild();
        // cout << "Number of remaining leafs = " <<
        // treeNew.pRoot->GetNumLeavesUnder() << endl; cout << "pRoot = " <<
        // (int) treeNew.pRoot << endl;
      }
    }

    pLeaf->RemoveLeafSelf();
    // cout << "After removing self\n";
    delete pLeaf;
    pLeaf = NULL;
    // cout << "here0\n";
  }
  // cout << "here\n";
  // cout << "Edge ids = ";
  // DumpIntVec( leavesEdgeIndices );
  // now we have the id we want as follows
  int res = 0;
  for (int lv = 2; lv < numLeaves; ++lv) {
    int base = 2 * lv - 1;
    res = res * base + leavesEdgeIndices[lv];
  }
  // cout << "res = " << res << endl;
  return res;
}

bool RBT ::RemoveLeaf(int lvid) {
  // first find the leaf
  int dummy;
  RBTNode *plf = pRoot->FindLeaf(lvid, dummy);
  if (plf == NULL) {
    // can not find the leaf
    return false;
  }

  // caution: if the leave is dirctly under root. then we have to change ROOT!
  if (plf->GetParent() == this->pRoot) {
    // set root to the sibling
    this->pRoot = plf->GetSibling();
    YW_ASSERT_INFO(this->pRoot != NULL, "Wrong: root becomes bad!");
  }

  plf->RemoveLeafSelf();
  // delete plf;
  plf = NULL;
  return true;
}

bool RBT ::IsSame(const RBT &tr) const {
  // test two trees are equivalent or not
  string trs = tr.GetNewick();
  string s0 = GetNewick();
  // when the leaf are ordered in a specific way,
  // two RBTs are the same iff Newick string is the same
  return trs == s0;
}
string RBT ::GetNewick() const {
  YW_ASSERT_INFO(pRoot != NULL, "Fail");
  return pRoot->GetNewick();
}

void RBT ::PruneLargeIdNodes(int idThres) {
  // get rid of id that is too large. possibly due to ARG issue
  // simply do an iteration
  TraversRecord tr;
  InitPostorderTranvers(tr);
  while (true) {
    //
    if (tr.pCurNode->IsLeaf() == true) {
      if (tr.pCurNode->GetLeafId() >= idThres) {
        // update current node
        RBTNode *pn = tr.pCurNode;
        RBTNode *pParNodeRem = pn->GetParent(); // the parent node is also gone
        // remove it
        NextPostorderTranvers(tr);
        if (tr.pCurNode == pParNodeRem) {
          NextPostorderTranvers(tr);
        }
        // cout << "Node extra removed: "  << pn->GetLeafId() << endl;
        pn->RemoveLeafSelf();
        // delete pn;
        pn = NULL;
        continue;
      }
    }

    // continue
    if (NextPostorderTranvers(tr) == false) {
      break;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////

bool RBT ::InitPostorderTranvers(TraversRecord &tr) {
  YW_ASSERT_INFO(pRoot != NULL, "Tree not initialized");

  // move down to the left-most leave (should be 0, verify it)
  RBTNode *pcur = this->pRoot->GetLeftMostChild();
  //    YW_ASSERT_INFO( pcur->GetLeafId() == 0, "The leftmost leaf must be 0" );
  tr.pCurNode = pcur;
  return true;
}

bool RBT ::NextPostorderTranvers(TraversRecord &tr) {
  // if we are at the root, we are done
  RBTNode *pCur = tr.pCurNode;
  if (pCur->GetParent() == NULL) {
    return false;
  }

  // if this is the left child, now move to right
  if (pCur->IsLeftChild() == true) {
    // start still from the left leaf
    tr.pCurNode = pCur->GetParent()->GetRightChild()->GetLeftMostChild();
  } else {
    // if it is right child, move up
    tr.pCurNode = pCur->GetParent();
  }
  return true;
}

void RBT ::RetrieveSplits() {
  // find and store all splits
  // we do this by retrieving splits in it
  // note we only store one side of splits, which contains 0
  TraversRecord tr;
  InitPostorderTranvers(tr);
  while (true) {
    set<int> lvs;
    tr.pCurNode->GetLeaves(lvs);
    if (lvs.find(0) != lvs.end()) {
      // save it
      if ((int)lvs.size() < this->numLeaves) {
        mapSplitsInTree.insert(map<set<int>, bool>::value_type(lvs, true));
      }
    } else {
      // store its complement
      set<int> compls;
      PopulateSetWithInterval(compls, 0, numLeaves - 1);
      SubtractSets(compls, lvs);
      if ((int)lvs.size() < this->numLeaves) {
        mapSplitsInTree.insert(map<set<int>, bool>::value_type(compls, true));
      }
    }

    // move to the next
    if (NextPostorderTranvers(tr) == false) {
      break;
    }
  }
}

RBTNode *RBT ::FindLeaf(int lvidParm, int &ponid) {
  // cout << "FindLeaf: lvidParm = " << lvidParm << endl;
  // just delegate to the root
  return this->pRoot->FindLeaf(lvidParm, ponid);
}

void RBT ::GetLeaves(set<int> &lvs) { pRoot->GetLeaves(lvs); }

void RBT ::Dump() const {
  pRoot->Dump();
  cout << endl;
}

void RBT ::DeleteLeaves(set<int> &lvids) {
  // delete leaves designated
  // here is a DUMB method: remove one by one
  // SLOW! but maybe enough for now. TBD
  for (set<int>::iterator it = lvids.begin(); it != lvids.end(); ++it) {
    int id = *it;
    // int dummy;
    // find the leave
    // RBTNode *tnode = FindLeaf( id, dummy);
    // if( tnode == NULL )
    //{
    //	cout << "Warning: leave id = " << id << " is not in the tree.\n";
    //	continue;
    //}
    // remove it
    if (RemoveLeaf(id) == false) {
      cout << "Warning: leave id = " << id << " is not in the tree.\n";
    }
    // cout << "After deleting leave = " << id << ", tree becomes: ";
    // Dump();
  }
}

void RBT ::RealignLeaves() {
  // cout << "RealignLeaves: tree = ";
  // Dump();
  // sometimes, say after leave is deleted, leaves are no longer contiguous,
  // this op sets it back to contiguous get all livids first
  set<int> lvids;
  GetLeaves(lvids);
  // convert to a lookup map
  map<int, int> mapLvidToRank;
  int rank = 0;
  for (set<int>::iterator it = lvids.begin(); it != lvids.end(); ++it) {
    mapLvidToRank.insert(map<int, int>::value_type(*it, rank++));
  }

  // now traversal the tree and do a traversal and reset leaf ids
  TraversRecord tr1;
  InitPostorderTranvers(tr1);
  while (true) {
    if (tr1.pCurNode->IsLeaf() == true) {
      int id = tr1.pCurNode->GetLeafId();
      YW_ASSERT_INFO(mapLvidToRank.find(id) != mapLvidToRank.end(),
                     "Leaf must be present");
      tr1.pCurNode->SetLeafId(mapLvidToRank[id]);
    }
    if (NextPostorderTranvers(tr1) == false) {
      break;
    }
  }

  // cout << "RealignLeaves: after realign, tree = ";
  // Dump();

  // also here also readjust the number of leaves
  this->numLeaves = lvids.size();
}

///////////////////////////////////////////////////////////////////////////////
// Other type of reconstruction

bool RBT ::ReconstructByPlainDesc(const vector<int> &listNodeLabels,
                                  const vector<int> &listParentNodePos,
                                  const vector<double> &listEdgeDist) {
  YW_ASSERT_INFO(listNodeLabels.size() >= 3, "Too small a tree");

  // first step is to get, for each tree node, what are the two children
  // this helps us in reconstructing the RBT tree. NOTE, we only deal with
  // non-leaves
  int numTNodes = listNodeLabels.size();
  vector<int> listNodeLeftChild, listNodeRightChild;
  for (int i = 0; i < numTNodes; ++i) {
    // start with -1 to indicate they are not set
    listNodeLeftChild.push_back(-1);
    listNodeRightChild.push_back(-1);
  }
  for (int i = 0; i < numTNodes; ++i) {
    int ppos = listParentNodePos[i];
    // cout << "i = " << i << ", posi = " << ppos << endl;
    if (ppos < 0) {
      // must reach the root
      break;
    }

    if (listNodeLeftChild[ppos] < 0) {
      // save it
      listNodeLeftChild[ppos] = i;
    } else if (listNodeRightChild[ppos] < 0) {
      listNodeRightChild[ppos] = i;
    } else {
      YW_ASSERT_INFO(
          false,
          "The tree is not binary. We can only handle binary for now.\n");
    }
  }

  // cout << "Here..\n";
  // first clear the old tree if any
  if (pRoot != NULL) {
    pRoot->Clear();
    delete pRoot;
    pRoot = NULL;
  }
  // cout << "Here...\n";
  // now do for every node
  vector<RBTNode *> listRBTNodes;
  for (int i = 0; i < numTNodes; ++i) {
    YW_ASSERT_INFO((listNodeLeftChild[i] >= 0 && listNodeRightChild[i] >= 0) ||
                       (listNodeLeftChild[i] < 0 && listNodeRightChild[i] < 0),
                   "WRONG");
    // cout << "Adding node-" << i << endl;
    // if it is leaf
    if (listNodeLeftChild[i] < 0) {
      // cout << "A leaf\n";
      //
      RBTNode *pn0 = new RBTNode(i);

      // also set height to be 1.0
      pn0->SetHeight(1.0);

      //
      listRBTNodes.push_back(pn0);
    } else {
      // cout << "Not a leaf\n";
      // not leaves
      int pnLeftInd = listNodeLeftChild[i];
      YW_ASSERT_INFO(pnLeftInd < numTNodes, "Tree node indices wrong");
      int pnRightInd = listNodeRightChild[i];
      YW_ASSERT_INFO(pnRightInd < numTNodes, "Tree node indices wrong");
      RBTNode *pn0 = listRBTNodes[pnLeftInd];
      RBTNode *pn1 = listRBTNodes[pnRightInd];
      RBTNode *prn;
      if (pn0->GetMinLeaveId() < pn1->GetMinLeaveId()) {
        prn = new RBTNode(pn0, pn1);
      } else {
        prn = new RBTNode(pn1, pn0);
      }

      // height is set by LEFT node (this may cause problem when the right is
      // NOT consistent for now, we IGNORE this hazard. NOTE, the leaf is
      // LOWEST, so we decrease
      double ht = pn0->GetHeight() - listEdgeDist[pnLeftInd];
      prn->SetHeight(ht);

      //
      listRBTNodes.push_back(prn);
    }
  }

  // set root
  int numNodesInList = listRBTNodes.size();
  YW_ASSERT_INFO(numNodesInList == numTNodes,
                 "Wrong in ReconstructByPlainDesc");
  this->pRoot = listRBTNodes[numNodesInList - 1];

  return true;
}

void RBT ::RetrievePlainDesc(int &numLvs, vector<int> &listNodeLabels,
                             vector<int> &listParentNodePos,
                             vector<double> &listEdgeDist) {
  numLvs = this->numLeaves;
  // init the return params
  listNodeLabels.clear();
  listParentNodePos.clear();
  listEdgeDist.clear();
  for (int i = 0; i < GetNodesNum(); ++i) {
    if (i < numLeaves) {
      listNodeLabels.push_back(i);
    } else {
      listNodeLabels.push_back(-1);
    }
    listParentNodePos.push_back(-1);
    listEdgeDist.push_back(-1.0);
  }

  // form a list of tree nodes
  // start iteration. Maintain TWO lists: one for leaves and one for internals
  vector<RBTNode *> listLeafNodes;
  // leaves are fixed
  listLeafNodes.resize(this->numLeaves);

  vector<RBTNode *> listInternalNodes;
  // use a map to quickly find location: well, not very elegent, but QUICK way
  // to do something
  map<RBTNode *, int> mapNodeToIndices;
  TraversRecord tr;
  InitPostorderTranvers(tr);
  while (true) {
    RBTNode *pcnode = tr.pCurNode;
    if (pcnode->IsLeaf() == true) {
      int lvid = pcnode->GetLeafId();
      YW_ASSERT_INFO(lvid >= 0 && lvid < this->numLeaves, "Fail in lvid");
      listLeafNodes[lvid] = pcnode;
      // save this node in the COMBINED index
      mapNodeToIndices.insert(map<RBTNode *, int>::value_type(pcnode, lvid));
    } else {
      listInternalNodes.push_back(pcnode);
      // save this node
      int ppos = listInternalNodes.size() - 1 + numLeaves;
      mapNodeToIndices.insert(map<RBTNode *, int>::value_type(pcnode, ppos));

      // also update input for the two children
      RBTNode *plc = pcnode->GetLeftChild();
      RBTNode *prc = pcnode->GetRightChild();
      YW_ASSERT_INFO(mapNodeToIndices.find(plc) != mapNodeToIndices.end(),
                     "WRONG");
      YW_ASSERT_INFO(mapNodeToIndices.find(prc) != mapNodeToIndices.end(),
                     "WRONG");
      int plcind = mapNodeToIndices[plc];
      int prcind = mapNodeToIndices[prc];
      // cout << "set left child node " << plcind << " to " << ppos << endl;
      // cout << "set right child node " << prcind << " to " << ppos << endl;
      listParentNodePos[plcind] = ppos;
      listParentNodePos[prcind] = ppos;

      // set edge length too
      double htPar = pcnode->GetHeight();
      double plcHt = plc->GetHeight();
      double prcHt = prc->GetHeight();
      if (htPar < 0 || plcHt < 0 || prcHt < 0 || plcHt < htPar ||
          prcHt < htPar) {
        // NOT VERY GOOD. TBD. 100707
        // set some arbitary number
        listEdgeDist[plcind] = 0.0;
        listEdgeDist[prcind] = 0.0;
      } else {
        // YW_ASSERT_INFO(htPar >= 0.0, "Height not set.");
        YW_ASSERT_INFO(plcHt >= 0.0 && prcHt >= 0.0, "Height not set.");
        YW_ASSERT_INFO(plcHt >= htPar && prcHt >= htPar, "Height not set.");
        listEdgeDist[plcind] = plcHt - htPar;
        listEdgeDist[prcind] = prcHt - htPar;
      }
    }

    //
    if (NextPostorderTranvers(tr) == false) {
      break;
    }
  }

  // make sure everything is correct
  for (int i = 0; i < GetNodesNum() - 1; ++i) {
    YW_ASSERT_INFO(listParentNodePos[i] >= 0 && listEdgeDist[i] >= 0,
                   "Some nodes are not correctly set.");
  }
}

void RBT ::AugamentDupRows(const vector<REMOVED_ROWS_INFO> &rmLvsStage) {
  // cout << "Before row augament, tree = ";
  // Dump();

  // ASSUMPTION: the tree is currently labeled from 0 - numLeaves-1
  // restore the leaves removed during matrix preprocessing
  // the list of items is taken out in a step-by-step procedure
  // IMPORTANT: need to reverse, since we start from the removed items
  for (int i = (int)rmLvsStage.size() - 1; i >= 0; --i) {
    // first reset the ids of the leaves
    // this is how the new ids will be
    int curLvNum = this->numLeaves;
    // int numRemoved = rmLvsStage[i].rowsRemoved.size();
    vector<int> vecRemRows;
    PopulateVecBySet(vecRemRows, rmLvsStage[i].rowsRemoved);
    vector<int> listOrigLeaveIds;
    GetOrigPositionAfterRemoval(curLvNum, vecRemRows, listOrigLeaveIds);
    // reconfig the leaves
    // cout << "Now setting leaves during tree augamentation...\n";
    SetLvids(listOrigLeaveIds);

    // then insert all the deleted rows back in
    // but first collect tips
    CollectTips();
    // now try to put back the removed rows
    for (int j = 0; j < (int)rmLvsStage[i].pairsRmKeepRows.size(); ++j) {
      //
      int rowNew = rmLvsStage[i].pairsRmKeepRows[j].first;
      YW_ASSERT_INFO(GetTip(rowNew) == NULL,
                     "Tip is already in"); // should not be already in
      int existId = rmLvsStage[i].pairsRmKeepRows[j].second;
      // cout << "existId = " << existId << ", rowNew = " << rowNew << endl;
      RBTNode *pn = GetTip(existId);
      YW_ASSERT_INFO(pn != NULL, "Src node not found");
      pn->AddSiblingToLeaf(rowNew);
      // cout << "After adding back " << rowNew << " the tree is: ";
      // Dump();
    }
    // update the number of leaves
    YW_ASSERT_INFO(rmLvsStage[i].pairsRmKeepRows.size() ==
                       rmLvsStage[i].rowsRemoved.size(),
                   "Removed record mismatch.");
    this->numLeaves += rmLvsStage[i].pairsRmKeepRows.size();
  }
  // cout << "After row augament, tree = ";
  // Dump();
}

void RBT ::SetLvids(const vector<int> &mapLvids) {
  // configure the name for the leaves
  // note, that we make assumption: the current leaves are labeled
  // consecutivatively!!!! OTHERWISE, it will not work wery well perform a
  // traversal
  TraversRecord tr;
  InitPostorderTranvers(tr);
  while (true) {
    if (tr.pCurNode->IsLeaf() == true) {
      // setup leave id
      int origId = tr.pCurNode->GetLeafId();
      YW_ASSERT_INFO(origId < (int)mapLvids.size(), "Leaf id is out of range");
      tr.pCurNode->SetLeafId(mapLvids[origId]);
      // cout << "Changing leave id from " << origId << " to " <<
      // mapLvids[origId] << endl;
    }
    if (NextPostorderTranvers(tr) == false) {
      break;
    }
  }
}

void RBT ::SetRoot(RBTNode *pRootNew) {
  // clear up if there is old root
  if (this->pRoot != NULL) {
    delete this->pRoot;
    this->pRoot = NULL;
  }
  YW_ASSERT_INFO(pRootNew != NULL, "Can not be NULL");
  this->pRoot = pRootNew;
  mapTipPtrs.clear();
  mapSplitsInTree.clear();
}
