#include "RerootTreeUtils.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <vector>

using namespace std;

void split(string &content, vector<string> &elements) {
  elements.clear();
  string tmp;
  for (int i = 0; i < (int)content.size(); ++i) {
    if (content[i] == ',' || content[i] == ':' || content[i] == '(' ||
        content[i] == ')') {
      if (!tmp.empty())
        elements.push_back(tmp);
      char ch[2] = { content[i], 0 };
      elements.push_back(string(ch));
      tmp.clear();
      continue;
    } else if ((content[i] >= '0' && content[i] <= '9') ||
               (content[i] >= 'A' && content[i] <= 'Z') ||
               (content[i] >= 'a' && content[i] <= 'z') || content[i] == '.') {
      const char ch[2] = { content[i], 0 };
      tmp.append(ch);
    }
  }
}
struct Edge {
  int a;
  double weight;
  Edge(int a, double weight) {
    this->a = a;
    this->weight = weight;
  }
  bool operator<(const Edge &edge) const {
    if (a != edge.a)
      return a < edge.a;
    return weight < edge.weight;
  }
};

double stringToDouble(string &content) {
  double ret = 0;
  int i = 0;
  for (; i < content.size() && content[i] != '.'; ++i) {
    if (content[i] < '0' || content[i] > '9') {
      printf("input tree string is not right\n");
      exit(0);
    }
    ret = ret * 10 + content[i] - '0';
  }
  double x = 0;
  if (content[i] == '.') {
    for (int j = content.size() - 1; j > i; --j) {
      if (content[j] < '0' || content[j] > '9') {
        printf("input tree string is not right\n");
        exit(0);
      }
      x = x * 0.1 + content[j] - '0';
    }
  }
  x = x * 0.1;
  return ret + x;
}

int stringToInt(string &content) {
  int ret = 0;
  for (int i = 0; i < content.size(); ++i)
    ret = ret * 10 + content[i] - '0';
  return ret;
}
void buildGraph(vector<string> &elements, map<int, map<int, double> > &graph,
                map<string, int> &leaf_to_label) {
  graph.clear();
  stack<char> s1;
  stack<Edge> s2;
  int a = -1;
  int cc = 0;
  for (int i = 0; i < elements.size(); ++i) {
    if (elements[i].compare("(") == 0) {
      s1.push('(');
    } else if (elements[i].compare(",") == 0) {
      s1.push(',');
    } else if (elements[i].compare(":") == 0) {
      s1.push(':');
    } else if (elements[i].compare(")") == 0) {
      if (s1.empty() || s1.top() != ',') {
        printf("input tree string is not right\n");
        exit(0);
      }
      s1.pop();
      if (s1.empty() || s1.top() != '(') {
        printf("input tree string is not right\n");
        exit(0);
      }
      s1.pop();
      a = cc;
      if ((int)s2.size() - 2 < 0) {
        printf("input tree string is not right\n");
        exit(0);
      }
      graph[a][s2.top().a] = s2.top().weight;
      graph[s2.top().a][a] = s2.top().weight;
      s2.pop();
      graph[a][s2.top().a] = s2.top().weight;
      graph[s2.top().a][a] = s2.top().weight;
      s2.pop();
      cc++;
    } else {
      if (s1.top() != ':') {
        a = cc;
        leaf_to_label[elements[i]] = cc;
        cc++;
      } else {
        double xx = stringToDouble(elements[i]);
        if (a == -1) {
          printf("input tree string is not right\n");
          exit(0);
        }
        s1.pop();
        s2.push(Edge(a, xx));
        a = -1;
      }
    }
  }
  if (!s1.empty() || !s2.empty()) {
    printf("input tree string is not right\n");
    exit(0);
  }
}

string convert(char *content, char *new_root) {
  string strRes;
  if (content == NULL || new_root == NULL)
    return strRes;
  string tree_str(content);
  vector<string> elements;
  split(tree_str, elements);
  map<int, map<int, double> > graph;
  map<string, int> leaf_to_label;
  buildGraph(elements, graph, leaf_to_label);

  string new_root_str(new_root);
  if (leaf_to_label.find(new_root_str) == leaf_to_label.end()) {
    printf("No such root %s\n", new_root);
    exit(0);
  }
  int nr = graph.size();
  int xx = -1;
  double yy = 0;

  // modify graph, add new root
  int nl = leaf_to_label[new_root_str];
  for (map<int, map<int, double> >::iterator iter = graph.begin();
       iter != graph.end(); ++iter) {
    if (iter->second.find(nl) != iter->second.end()) {
      yy = iter->second[nl];
      xx = iter->first;
      iter->second[nr] = yy / 2;
      iter->second.erase(nl);
      break;
    }
  }
  graph[nr][xx] = yy / 2;
  graph[nr][nl] = yy / 2;
  graph[nl].clear();
  graph[nl][nr] = yy / 2;
#if 0
printf("graph\n");
for (map<int, map<int, double> >::iterator iter1 =graph.begin();iter1!=graph.end();++iter1) {
printf("%d:", iter1->first);
for (map<int, double>::iterator iter2 =iter1->second.begin();iter2!=iter1->second.end();++iter2) {
printf("(%d,%lf) ", iter2->first, iter2->second);
}
printf("\n");
}
#endif
  // bfs, get new weight
  int n = graph.size();
  vector<double> wei;
  vector<bool> flag;
  wei.reserve(n);
  flag.reserve(n);
  for (int i = 0; i < n; i++) {
    flag.push_back(false);
    wei.push_back(0);
  }
  queue<int> qu;
  qu.push(nr);
  flag[nr] = true;
  map<int, set<int> > tree;
  map<int, int> parent;
  while (!qu.empty()) {
    int t = qu.front();
    qu.pop();
    if (graph.find(t) == graph.end())
      continue;
    for (map<int, double>::iterator iter = graph[t].begin();
         iter != graph[t].end(); ++iter) {
      if (flag[iter->first])
        continue;
      flag[iter->first] = true;
      qu.push(iter->first);
      wei[iter->first] = wei[t] + (iter->second);
      tree[t].insert(iter->first);
      parent[iter->first] = t;
    }
  }
#if 0
printf("tree\n");
for (map<int, set<int> >::iterator iter1 =tree.begin();iter1!=tree.end();++iter1) {
printf("%d:", iter1->first);
for (set<int>::iterator iter2 =iter1->second.begin();iter2!=iter1->second.end();++iter2) {
printf("%d ", *iter2);
}
printf("\n");
}
printf("parent\n");
for (map<int, int>::iterator iter =parent.begin();iter!=parent.end();++iter) {
printf("%d %d\n", iter->first, iter->second);
}


printf("weight\n");
for (int i=0;i<wei.size();++i) {
printf("%d %lf\n", i, wei[i]);
}
#endif
  // eliminate old root
  int old = n - 2;
  xx = *(tree[old].begin());
  for (map<int, set<int> >::iterator iter = tree.begin(); iter != tree.end();
       ++iter) {
    if (iter->second.find(old) != iter->second.end()) {
      iter->second.erase(old);
      iter->second.insert(xx);
      parent[xx] = iter->first;
      break;
    }
  }
  tree.erase(old);
  parent.erase(old);

#if 0
printf("tree\n");
for (map<int, set<int> >::iterator iter1 =tree.begin();iter1!=tree.end();++iter1) {
printf("%d:", iter1->first);
for (set<int>::iterator iter2 =iter1->second.begin();iter2!=iter1->second.end();++iter2) {
printf("%d ", *iter2);
}
printf("\n");
}
printf("parent\n");
for (map<int, int>::iterator iter =parent.begin();iter!=parent.end();++iter) {
printf("%d %d\n", iter->first, iter->second);
}


printf("weight\n");
for (int i=0;i<wei.size();++i) {
printf("%d %lf\n", i, wei[i]);
}
#endif
  // print new tree
  map<int, string> nts;

  deque<pair<int, int> > de;
  for (map<int, set<int> >::iterator iter = tree.begin(); iter != tree.end();
       ++iter) {
    if (iter->second.size() == 2) {
      int a[3] = { 0, 0, 0 };
      for (set<int>::iterator iter2 = iter->second.begin();
           iter2 != iter->second.end(); ++iter2) {
        if (tree.find(*iter2) == tree.end())
          a[++a[0]] = *iter2;
      }
      if (a[0] == 2) {
        de.push_back(pair<int, int>(a[1], a[2]));
      }
    }
  }

  for (map<string, int>::iterator iter = leaf_to_label.begin();
       iter != leaf_to_label.end(); ++iter) {
    char tmp[100];
    double tt = 0;
    if (parent.find(iter->second) != parent.end())
      tt = wei[parent[iter->second]];
    sprintf(tmp, "%f", wei[iter->second] - tt);
    if (iter->second != nr)
      nts[iter->second] = iter->first + ':' + string(tmp);
    else
      nts[iter->second] = iter->first;
  }
#if 0
printf("node to string\n");
for (map<int, string>::iterator iter =nts.begin();iter!=nts.end();++iter) {
printf("%d: %s\n", iter->first, iter->second.c_str());
}
#endif
  while (!de.empty()) {
    pair<int, int> a = de.front();
    de.pop_front();
    int pa = parent[a.first];
    char tmp[100];
    double tt = 0;
    if (parent.find(pa) != parent.end())
      tt = wei[parent[pa]];
    sprintf(tmp, "%f", wei[pa] - tt);
    if (pa != nr)
      nts[pa] =
          '(' + nts[a.first] + ',' + nts[a.second] + ')' + ':' + string(tmp);
    else
      nts[pa] = '(' + nts[a.first] + ',' + nts[a.second] + ')';
    tree.erase(pa);
    if (parent.find(pa) != parent.end()) {
      int ppa = parent[pa];
      int sibling;
      for (set<int>::iterator iter = tree[ppa].begin(); iter != tree[ppa].end();
           ++iter) {
        if ((*iter) != pa)
          sibling = *iter;
      }
      if (tree.find(sibling) == tree.end()) {
        de.push_back(pair<int, int>(pa, sibling));
      }
    }
  }
  // printf("%s\n", nts[nr].c_str());
  strRes = nts[nr];
  return strRes;
}

void Test_split() {
  string a("(((1:1.0,2:2.0):1.2,(3:1.0,4:2.0):1.6):1.5,5:1.0)");
  string b("( ( ( 1 : 1.0 , 2 : 2.0 ) : 1.2 , ( 3 : 1.0 , 4 : 2.0 ) : 1.6 ) : "
           "1.5 , 5 : 1.0 ) ");
  vector<string> elements;
  split(a, elements);
  for (int i = 0; i < elements.size(); ++i) {
    printf("%s ", elements[i].c_str());
  }
  printf("\n");
  split(b, elements);
  for (int i = 0; i < elements.size(); ++i) {
    printf("%s ", elements[i].c_str());
  }
  printf("\n");
}

void Test_buildGraph() {
  string a("(((1:1.0,2:2.0):1.2,(3:1.0,4:2.0):1.6):1.5,5:1.0)");
  vector<string> elements;
  split(a, elements);
  map<int, map<int, double> > graph;
  map<string, int> leaf_to_label;
  buildGraph(elements, graph, leaf_to_label);
  printf("leaf to label\n");
  for (map<string, int>::iterator iter = leaf_to_label.begin();
       iter != leaf_to_label.end(); ++iter) {
    printf("%s:%d\n", iter->first.c_str(), iter->second);
  }
  printf("Graph\n");
  for (map<int, map<int, double> >::iterator iter1 = graph.begin();
       iter1 != graph.end(); ++iter1) {
    printf("%d:", iter1->first);
    for (map<int, double>::iterator iter2 = iter1->second.begin();
         iter2 != iter1->second.end(); ++iter2) {
      printf("(%d,%lf) ", iter2->first, iter2->second);
    }
    printf("\n");
  }
}

string ReRootTreeNewick(char *nwFile, char *taxaNewRoot) {
  // char * a="(((1:1.0,2:2.0):1.2,(3:1.0,4:2.0):1.6):1.5,5:1.0)";
  // char *b ="3";
  // usage for converting
  return convert(nwFile, taxaNewRoot);
}
