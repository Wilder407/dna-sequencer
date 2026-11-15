// Needleman-Wunsch global sequence alignment
// Author: Lisa Wilder (GitHub: Wilder407)

#include "Matrix.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

static bool is_header(const string& line) {
  size_t k = line.find_first_not_of(" \t\r");
  if (k == string::npos) return false;
  return line[k] == '>' || line.compare(k, 2, "//") == 0;
}

string clean_sequence(const string& raw) {
  string out;
  istringstream in(raw);
  string line;
  while (getline(in, line)) {
    if (is_header(line)) continue;
    for (char c : line) {
      if (isspace(static_cast<unsigned char>(c)) || isdigit(static_cast<unsigned char>(c))) continue;
      out += static_cast<char>(toupper(static_cast<unsigned char>(c)));
    }
  }
  return out;
}

vector<string> read_sequences(const string& path) {
  ifstream f(path);
  if (!f) throw runtime_error("Could not open file: " + path);

  vector<string> seqs;
  string line, current;
  bool saw_header = false;
  while (getline(f, line)) {
    if (is_header(line)) {
      if (saw_header && !current.empty()) seqs.push_back(clean_sequence(current));
      current.clear();
      saw_header = true;
    } else {
      current += line + "\n";
    }
  }
  string last = clean_sequence(current);
  if (!last.empty()) seqs.push_back(last);
  return seqs;
}

bool valid_input(const string& seq) {
  if (seq.empty()) return false;
  return all_of(seq.begin(), seq.end(),
                [](char c) { return string("ACGTU").find(c) != string::npos; });
}

int match_score(char a, char b, const Scoring& s) {
  return a == b ? s.match : s.mismatch;
}

NWMatrices N_W(const string& x, const string& y, const Scoring& s) {
  const int m = x.size(), n = y.size();
  NWMatrices M;
  M.score.assign(m + 1, vector<int>(n + 1, 0));
  M.dir.assign(m + 1, vector<Dir>(n + 1, Dir::None));

  // First column/row: aligning i letters against nothing costs i gaps.
  // With free end gaps, leading gaps cost nothing, so these stay 0.
  const int edge = s.free_end_gaps ? 0 : s.gap;
  for (int i = 1; i <= m; i++) { M.score[i][0] = i * edge; M.dir[i][0] = Dir::Up; }
  for (int j = 1; j <= n; j++) { M.score[0][j] = j * edge; M.dir[0][j] = Dir::Left; }

  for (int i = 1; i <= m; i++) {
    for (int j = 1; j <= n; j++) {
      int diag = M.score[i - 1][j - 1] + match_score(x[i - 1], y[j - 1], s);
      int up   = M.score[i - 1][j] + s.gap;   // x[i-1] aligned to a gap
      int left = M.score[i][j - 1] + s.gap;   // y[j-1] aligned to a gap
      int best = max({diag, up, left});
      M.score[i][j] = best;
      // Tie-break order: diagonal, then up, then left.
      M.dir[i][j] = (best == diag) ? Dir::Diag : (best == up) ? Dir::Up : Dir::Left;
    }
  }
  return M;
}

int nw_score_only(const string& x, const string& y, const Scoring& s) {
  const int m = x.size(), n = y.size();
  const int edge = s.free_end_gaps ? 0 : s.gap;
  vector<int> prev(n + 1), cur(n + 1);
  int best_last_col = 0;  // best score in column n (free end gaps)
  for (int j = 0; j <= n; j++) prev[j] = j * edge;
  for (int i = 1; i <= m; i++) {
    cur[0] = i * edge;
    for (int j = 1; j <= n; j++) {
      cur[j] = max({prev[j - 1] + match_score(x[i - 1], y[j - 1], s),
                    prev[j] + s.gap, cur[j - 1] + s.gap});
    }
    if (i == 1 || cur[n] > best_last_col) best_last_col = cur[n];
    swap(prev, cur);
  }
  if (!s.free_end_gaps) return prev[n];
  int best = *max_element(prev.begin(), prev.end());  // last row
  return (m > 0) ? max(best, best_last_col) : best;
}

Alignment traceback(const NWMatrices& M, const string& x, const string& y,
                    const Scoring& s) {
  Alignment a;
  const int m = x.size(), n = y.size();
  int i = m, j = n;

  if (s.free_end_gaps) {
    // Start from the best cell in the last row or column; whatever is left
    // of the other sequence becomes free trailing gaps.
    int bi = m, bj = n;
    for (int jj = 0; jj <= n; jj++) if (M.score[m][jj] > M.score[bi][bj]) { bi = m; bj = jj; }
    for (int ii = 0; ii <= m; ii++) if (M.score[ii][n] > M.score[bi][bj]) { bi = ii; bj = n; }
    for (int jj = n; jj > bj; jj--) { a.aligned_x += '-'; a.aligned_y += y[jj - 1]; a.gaps++; }
    for (int ii = m; ii > bi; ii--) { a.aligned_x += x[ii - 1]; a.aligned_y += '-'; a.gaps++; }
    i = bi; j = bj;
  }
  a.score = M.score[i][j];

  while (i > 0 || j > 0) {
    Dir d = M.dir[i][j];
    if (d == Dir::Diag) {
      a.aligned_x += x[i - 1];
      a.aligned_y += y[j - 1];
      (x[i - 1] == y[j - 1]) ? a.matches++ : a.mismatches++;
      i--; j--;
    } else if (d == Dir::Up) {
      a.aligned_x += x[i - 1];
      a.aligned_y += '-';
      a.gaps++;
      i--;
    } else {
      a.aligned_x += '-';
      a.aligned_y += y[j - 1];
      a.gaps++;
      j--;
    }
  }
  // Built back-to-front; reverse once instead of prepending each step.
  reverse(a.aligned_x.begin(), a.aligned_x.end());
  reverse(a.aligned_y.begin(), a.aligned_y.end());
  a.length = a.aligned_x.size();

  // End gaps: gap columns before the first and after the last paired column.
  auto paired = [&](int k) { return a.aligned_x[k] != '-' && a.aligned_y[k] != '-'; };
  int first = 0, last = a.length - 1;
  while (first < a.length && !paired(first)) first++;
  while (last >= 0 && !paired(last)) last--;
  a.end_gaps = (first >= a.length) ? a.length : first + (a.length - 1 - last);
  return a;
}

static const char* arrow(Dir d) {
  switch (d) {
    case Dir::Diag: return "↖";
    case Dir::Up:   return "↑";
    case Dir::Left: return "←";
    default:        return "·";
  }
}

void print_matrix(const NWMatrices& M, const string& x, const string& y) {
  auto header = [&]() {
    cout << "\t\t-";
    for (char c : y) cout << "\t" << c;
    cout << "\n";
  };

  cout << "Score matrix (F):\n";
  header();
  for (size_t i = 0; i < M.score.size(); i++) {
    cout << "\t" << (i == 0 ? '-' : x[i - 1]);
    for (int v : M.score[i]) cout << "\t" << v;
    cout << "\n";
  }

  cout << "\nDirection matrix:\n";
  header();
  for (size_t i = 0; i < M.dir.size(); i++) {
    cout << "\t" << (i == 0 ? '-' : x[i - 1]);
    for (Dir d : M.dir[i]) cout << "\t" << arrow(d);
    cout << "\n";
  }
}

void print_alignment(const Alignment& a, int width) {
  string mid;
  for (int k = 0; k < a.length; k++) {
    char p = a.aligned_x[k], q = a.aligned_y[k];
    mid += (p == '-' || q == '-') ? ' ' : (p == q ? '|' : '.');
  }
  for (int k = 0; k < a.length; k += width) {
    cout << "x  " << a.aligned_x.substr(k, width) << "\n"
         << "   " << mid.substr(k, width) << "\n"
         << "y  " << a.aligned_y.substr(k, width) << "\n\n";
  }
}
