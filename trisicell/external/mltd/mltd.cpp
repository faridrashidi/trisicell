#include "mltd.h"

using std::cerr;
using std::cout;
using std::ifstream;
using std::make_pair;
using std::map;
using std::max;
using std::pair;
using std::string;
using std::stringstream;
using std::vector;

struct Edge {
  int s, t, c, f, a;
};

struct Graph {
  int cedge = 0;
  vector<Edge> edges;
  int n;
  vector<vector<int> > g;
  void init(int _n) {
    n = _n;
    g = vector<vector<int> >(n, vector<int>());
  }
  void addEdge(int u, int v, int c, int a) {
    edges.push_back(Edge{ u, v, c, 0, a });
    g[u].push_back(cedge++);
    edges.push_back(Edge{ v, u, 0, 0, -a });
    g[v].push_back(cedge++);
  }
};

map<string, int> label2int;
int nlabel = 0;

struct Tree {
  vector<vector<int> > g;
  vector<vector<int> > l;
  map<pair<int, int>, int> pair2int;
  int npair = 0;
  int nname = 0;
  int r = -1;
  int nl = 0;
  vector<int> p;
  void dfs(int p, int s) {
    pair2int[make_pair(p, s)] = npair++;
    for (auto x : g[s]) {
      dfs(p, x);
    }
  }
  void read(const char *file) {
    map<string, int> name2int;
    ifstream in(file);
    string str;
    while (getline(in, str)) {
      stringstream ss;
      int flag = -1;
      for (auto &c : str) {
        if (c == ':') {
          flag = 0;
          ss << ' ';
        } else if (c == '=') {
          flag = 1;
          ss << ' ';
        } else if (c == ',') {
          ss << ' ';
        } else {
          ss << c;
        }
      }
      string s;
      assert(flag == 1 || flag == 0);
      if (flag == 1) {
        ss >> s;
        if (!name2int.count(s)) {
          name2int[s] = nname++;
          l.push_back(vector<int>());
          p.push_back(-1);
        }
        int x = name2int[s];
        while (ss >> s) {
          if (!label2int.count(s)) {
            label2int[s] = nlabel++;
          }
          l[x].push_back(label2int[s]);
          nl++;
        }
      }
      if (flag == 0) {
        ss >> s;
        if (!name2int.count(s)) {
          name2int[s] = nname++;
          l.push_back(vector<int>());
          p.push_back(-1);
        }
        int x = name2int[s];
        while (ss >> s) {
          if (!name2int.count(s)) {
            name2int[s] = nname++;
            l.push_back(vector<int>());
            p.push_back(-1);
          }
          p[name2int[s]] = x;
        }
      }
    }
    int cnt = 0;
    g.resize(p.size());
    // cerr << p.size() << "\n";
    for (int i = 0; i < (int)p.size(); i++) {
      cnt += p[i] == -1;
      if (p[i] != -1) {
        g[p[i]].push_back(i);
      } else {
        r = i;
      }
    }
    assert(cnt == 1);
    for (int i = 0; i < (int)g.size(); i++) {
      dfs(i, i);
    }
  }
};

vector<vector<int> > T;
vector<vector<int> > D;
vector<vector<int> > G;
Tree t, s;

int cD(int a, int b, int c, int d) {
  int &res = D[t.pair2int[make_pair(a, c)]][s.pair2int[make_pair(b, d)]];
  if (res != -1)
    return res;
  if (a == c && b == d) {
    return res = T[c][d];
  }
  if (a == c && b != d) {
    return res = T[c][d] + cD(a, b, c, s.p[d]);
  }
  if (a != c && b == d) {
    return res = T[c][d] + cD(a, b, t.p[c], d);
  }
  return res = T[c][d] + max(cD(a, b, t.p[c], d), cD(a, b, c, s.p[d]));
}

vector<int> subtree(Tree &g, int s) {
  vector<int> res;
  res.push_back(s);
  for (auto x : g.g[s]) {
    auto q = subtree(g, x);
    for (auto y : q) {
      res.push_back(y);
    }
  }
  return res;
}

int mcmf(Graph &g, int S, int T) {
  vector<int> d(g.n);
  vector<int> p(g.n);
  int res = 0;
  int inf = std::numeric_limits<int>::max() / 2;
  do {
    std::fill(d.begin(), d.end(), inf);
    std::fill(p.begin(), p.end(), -1);
    d[S] = 0;
    int found;
    do {
      found = 0;
      for (int i = 0; i < g.n; i++) {
        if (d[i] == inf)
          continue;
        for (auto id : g.g[i]) {
          auto s = g.edges[id].s;
          auto t = g.edges[id].t;
          auto c = g.edges[id].c;
          auto f = g.edges[id].f;
          auto a = g.edges[id].a;
          if (c - f > 0 && d[s] + a < d[t]) {
            found = 1;
            d[t] = d[s] + a;
            p[t] = id;
          }
        }
      }
    } while (found);
    if (d[T] > 0)
      break;
    int v = T;
    res += d[T];
    while (p[v] != -1) {
      auto id = p[v];
      g.edges[id].f++;
      g.edges[id ^ 1].f--;
      v = g.edges[id].s;
    }
  } while (1);
  return res;
}

int cG(int a, int b) {
  int &res = G[a][b];
  if (res != -1)
    return res;
  int na = t.g[a].size();
  int nb = s.g[b].size();
  Graph g;
  g.init(na + nb + 2);
  int S = na + nb;
  int T = S + 1;
  for (int i = 0; i < na; i++) {
    g.addEdge(S, i, 1, 0);
  }
  for (int i = 0; i < nb; i++) {
    g.addEdge(na + i, T, 1, 0);
  }
  auto ca = vector<vector<int> >(t.g[a].size());
  auto cb = vector<vector<int> >(s.g[b].size());
  for (int i = 0; i < (int)t.g[a].size(); i++) {
    ca[i] = subtree(t, t.g[a][i]);
  }
  for (int i = 0; i < (int)s.g[b].size(); i++) {
    cb[i] = subtree(s, s.g[b][i]);
  }
  for (int i = 0; i < (int)ca.size(); i++) {
    for (int j = 0; j < (int)cb.size(); j++) {
      int t = 0;
      for (auto x : ca[i]) {
        for (auto y : cb[j]) {
          t = max(t, cD(ca[i][0], cb[j][0], x, y) + cG(x, y));
        }
      }
      g.addEdge(i, na + j, 1, -t);
    }
  }
  return res = -mcmf(g, S, T);
}

MLTDResult calc_mltd(const char *tree1, const char *tree2) {
  t.read(tree1);
  s.read(tree2);
  auto tl = vector<int>(nlabel, -1);
  auto sl = vector<int>(nlabel, -1);
  for (int i = 0; i < (int)t.l.size(); i++) {
    for (auto &x : t.l[i]) {
      assert(tl[x] == -1);
      tl[x] = i;
    }
  }
  for (int i = 0; i < (int)s.l.size(); i++) {
    for (auto &x : s.l[i]) {
      assert(sl[x] == -1);
      sl[x] = i;
    }
  }
  T = vector<vector<int> >(t.nname, vector<int>(s.nname, 0));
  for (int i = 0; i < nlabel; i++) {
    if (sl[i] != -1 && tl[i] != -1) {
      T[tl[i]][sl[i]]++;
    }
  }
  /*for (auto &x : T) {
    for (auto &y : x) {
      cerr << y << " ";
    }
    cerr << "\n";
  }
  cerr << t.npair << " " << s.npair << "\n";*/
  D = vector<vector<int> >(t.npair, vector<int>(s.npair, -1));
  G = vector<vector<int> >(t.nname, vector<int>(s.nname, -1));
  int res = 0;
  // cerr << t.r << " " << s.r << "\n";
  for (int i = 0; i < t.nname; i++) {
    for (int j = 0; j < s.nname; j++) {
      res = max(res, cG(i, j) + cD(t.r, s.r, i, j));
    }
  }
  /*cerr << "T1\n";
  for (int i = 0; i < t.nname; i++) {
    cerr << i << " :(";
    for (auto x : t.l[i]) {
      cerr << x << ",";
    }
    cerr << ") ";
    for (auto x : t.g[i]) {
      cerr << x << " ";
    }
    cerr << "\n";
  }
  cerr << "T2\n";
  for (int i = 0; i < s.nname; i++) {
    cerr << i << " :(";
    for (auto x : s.l[i]) {
      cerr << x << ",";
    }
    cerr << ") ";
    for (auto x : s.g[i]) {
      cerr << x << " ";
    }
    cerr << "\n";
  }*/
  MLTDResult mltdresult;
  // cout << "\nOutput:\n\n";
  // cout << "Distance = " << t.nl + s.nl - 2 * res << "\n";
  // cout << "Similarity = " << res << "\n";
  // cout << "Normalized Similarity = " << res * 1.0 / nlabel << "\n\n";
  mltdresult.distance = t.nl + s.nl - 2 * res;
  mltdresult.similarity = res;
  mltdresult.normalized_similarity = res * 1.0 / nlabel;
  return mltdresult;
}

int main(int argc, char *argv[]) { return 0; }
