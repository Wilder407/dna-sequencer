// Unit tests for the Needleman-Wunsch implementation.
// Build and run with: make test

#include "Matrix.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace std;

static int failures = 0, checks = 0;

#define CHECK(cond)                                                     \
  do {                                                                  \
    checks++;                                                           \
    if (!(cond)) {                                                      \
      failures++;                                                       \
      cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  " #cond "\n"; \
    }                                                                   \
  } while (0)

// Removes '-' so an aligned row can be compared with its input sequence.
static string ungap(const string& s) {
  string out;
  for (char c : s) if (c != '-') out += c;
  return out;
}

// Rescores an alignment column by column.
static int rescore(const Alignment& a, const Scoring& s) {
  int total = 0;
  for (int k = 0; k < a.length; k++) {
    char p = a.aligned_x[k], q = a.aligned_y[k];
    total += (p == '-' || q == '-') ? s.gap : match_score(p, q, s);
  }
  return total;
}

int main() {
  // Textbook example (Wikipedia, "Needleman-Wunsch algorithm"):
  // match +1, mismatch -1, gap -1 gives an optimal score of 0.
  {
    Scoring s{1, -1, -1};
    auto M = N_W("GATTACA", "GCATGCU", s);
    Alignment a = traceback(M, "GATTACA", "GCATGCU", s);
    CHECK(a.score == 0);
    CHECK(rescore(a, s) == a.score);
    CHECK(ungap(a.aligned_x) == "GATTACA");
    CHECK(ungap(a.aligned_y) == "GCATGCU");
  }

  // Identical sequences: every column is a match.
  {
    Scoring s;
    auto M = N_W("ACGTACGT", "ACGTACGT", s);
    Alignment a = traceback(M, "ACGTACGT", "ACGTACGT", s);
    CHECK(a.score == 8);
    CHECK(a.matches == 8 && a.mismatches == 0 && a.gaps == 0);
    CHECK(a.identity() == 100.0);
  }

  // One deletion: ACGT vs AGT should align as ACGT / A-GT.
  {
    Scoring s;  // +1 / -1 / -2
    auto M = N_W("ACGT", "AGT", s);
    Alignment a = traceback(M, "ACGT", "AGT", s);
    CHECK(a.score == 1);
    CHECK(a.aligned_x == "ACGT");
    CHECK(a.aligned_y == "A-GT");
    CHECK(a.gaps == 1);
  }

  // First row and column hold multiples of the gap penalty.
  {
    Scoring s;
    auto M = N_W("AC", "GTA", s);
    CHECK(M.score[0][0] == 0);
    CHECK(M.score[2][0] == 2 * s.gap);
    CHECK(M.score[0][3] == 3 * s.gap);
  }

  // Score-only version agrees with the full matrix, and the traceback
  // always rescores to F[m][n], on random sequences.
  {
    mt19937 rng(7);
    const string bases = "ACGT";
    for (int t = 0; t < 200; t++) {
      string x, y;
      int lx = rng() % 30 + 1, ly = rng() % 30 + 1;
      for (int k = 0; k < lx; k++) x += bases[rng() % 4];
      for (int k = 0; k < ly; k++) y += bases[rng() % 4];
      Scoring s{2, -1, -3};
      auto M = N_W(x, y, s);
      Alignment a = traceback(M, x, y, s);
      CHECK(nw_score_only(x, y, s) == a.score);
      CHECK(rescore(a, s) == a.score);
      CHECK(ungap(a.aligned_x) == x && ungap(a.aligned_y) == y);
    }
  }

  // Free end gaps: a partial read inside a longer sequence scores as if
  // the overhangs were not there.
  {
    Scoring s;
    s.free_end_gaps = true;
    const string full = "TTTTACGTACGTGGGG", part = "ACGTACGT";
    auto M = N_W(full, part, s);
    Alignment a = traceback(M, full, part, s);
    CHECK(a.score == 8);
    CHECK(a.matches == 8 && a.mismatches == 0);
    CHECK(a.end_gaps == 8 && a.gaps == 8);
    CHECK(nw_score_only(full, part, s) == 8);
    CHECK(ungap(a.aligned_x) == full && ungap(a.aligned_y) == part);

    Scoring g;  // global mode charges for the same overhangs
    CHECK(nw_score_only(full, part, g) == 8 + 8 * g.gap);
  }

  // Score-only and full matrix agree in free-end-gap mode on random input.
  {
    mt19937 rng(11);
    const string bases = "ACGT";
    for (int t = 0; t < 200; t++) {
      string x, y;
      int lx = rng() % 30 + 1, ly = rng() % 30 + 1;
      for (int k = 0; k < lx; k++) x += bases[rng() % 4];
      for (int k = 0; k < ly; k++) y += bases[rng() % 4];
      Scoring s{1, -1, -2, true};
      auto M = N_W(x, y, s);
      Alignment a = traceback(M, x, y, s);
      CHECK(nw_score_only(x, y, s) == a.score);
      CHECK(ungap(a.aligned_x) == x && ungap(a.aligned_y) == y);
    }
  }

  // Input cleaning and validation.
  {
    CHECK(clean_sequence(">seq1\nacg t\n  1 gca\n") == "ACGTGCA");
    CHECK(valid_input("ACGTU"));
    CHECK(!valid_input(""));
    CHECK(!valid_input("ACGN"));
  }

  // Reading two records from one file, with '//' headers like the data file.
  {
    const char* path = "test_tmp_seqs.txt";
    ofstream(path) << "\n// first\nacgt acgt\n\n// second\nggcc\n";
    auto seqs = read_sequences(path);
    remove(path);
    CHECK(seqs.size() == 2);
    CHECK(seqs.size() == 2 && seqs[0] == "ACGTACGT" && seqs[1] == "GGCC");
  }

  cout << (checks - failures) << "/" << checks << " checks passed\n";
  return failures ? 1 : 0;
}
