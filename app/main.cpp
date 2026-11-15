// Needleman-Wunsch global sequence alignment: command-line program
// Author: Lisa Wilder (GitHub: Wilder407)
//
// The original version (CSCI 2270, spring 2025) was written with help from
// ChatGPT-4o (OpenAI), April 27 and May 4, 2025.

#include "Matrix.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace std;

static void usage(const char* prog) {
  cerr << "Usage: " << prog << " [options] [FILE1 [FILE2]]\n\n"
       << "  FILE1 FILE2       one sequence per file (FASTA or plain text)\n"
       << "  FILE1             a file with two FASTA records\n"
       << "  (no files)        type or paste the two sequences\n\n"
       << "Options:\n"
       << "  --match N         score for a match (default 1)\n"
       << "  --mismatch N      score for a mismatch (default -1)\n"
       << "  --gap N           gap penalty per gap position (default -2)\n"
       << "  --free-end-gaps   don't penalize gaps at either end (semi-global);\n"
       << "                    use when one sequence is a partial read of the other\n"
       << "  --show-matrix     print the score and direction matrices\n"
       << "                    (printed automatically when both sequences are <= 15)\n"
       << "  --shuffles N      compare the score to N shuffled versions of y (default 0)\n"
       << "  --seed N          random seed for --shuffles (default 42)\n";
}

static string prompt_sequence(const string& label) {
  cout << "Enter the " << label << " sequence (end with a blank line):\n";
  string raw, line;
  while (getline(cin, line) && !line.empty()) raw += line + "\n";
  return clean_sequence(raw);
}

int main(int argc, char* argv[]) {
  Scoring s;
  bool show_matrix = false;
  int shuffles = 0;
  unsigned seed = 42;
  vector<string> files;

  for (int k = 1; k < argc; k++) {
    string arg = argv[k];
    auto next_int = [&](int& out) {
      if (k + 1 >= argc) { usage(argv[0]); exit(1); }
      out = stoi(argv[++k]);
    };
    if (arg == "--match") next_int(s.match);
    else if (arg == "--mismatch") next_int(s.mismatch);
    else if (arg == "--gap") next_int(s.gap);
    else if (arg == "--shuffles") next_int(shuffles);
    else if (arg == "--seed") { int v; next_int(v); seed = v; }
    else if (arg == "--show-matrix") show_matrix = true;
    else if (arg == "--free-end-gaps") s.free_end_gaps = true;
    else if (arg == "-h" || arg == "--help") { usage(argv[0]); return 0; }
    else if (!arg.empty() && arg[0] == '-') { usage(argv[0]); return 1; }
    else files.push_back(arg);
  }

  string x, y;
  try {
    if (files.size() == 2) {
      auto a = read_sequences(files[0]);
      auto b = read_sequences(files[1]);
      if (a.empty() || b.empty()) throw runtime_error("A file contained no sequence.");
      x = a[0];
      y = b[0];
    } else if (files.size() == 1) {
      auto seqs = read_sequences(files[0]);
      if (seqs.size() < 2) throw runtime_error("Expected two sequences in " + files[0]);
      x = seqs[0];
      y = seqs[1];
    } else if (files.empty()) {
      cout << "Needleman-Wunsch global alignment of two nucleotide sequences\n\n";
      x = prompt_sequence("first");
      y = prompt_sequence("second");
    } else {
      usage(argv[0]);
      return 1;
    }
  } catch (const exception& e) {
    cerr << "Error: " << e.what() << "\n";
    return 1;
  }

  if (!valid_input(x) || !valid_input(y)) {
    cerr << "Error: sequences must be non-empty and use only A, C, G, T, U.\n";
    return 1;
  }

  NWMatrices M = N_W(x, y, s);
  Alignment a = traceback(M, x, y, s);

  cout << "Scoring: match " << s.match << ", mismatch " << s.mismatch
       << ", gap " << s.gap
       << (s.free_end_gaps ? " (end gaps free)" : "") << "\n";
  cout << "Length x: " << x.size() << "   Length y: " << y.size() << "\n\n";

  if (show_matrix || (x.size() <= 15 && y.size() <= 15)) {
    print_matrix(M, x, y);
    cout << "\n";
  }

  print_alignment(a);

  cout << fixed << setprecision(1);
  cout << "Alignment score: " << a.score << "\n"
       << "Alignment length: " << a.length << " columns\n"
       << "Matches: " << a.matches << "   Mismatches: " << a.mismatches
       << "   Gap positions: " << a.gaps << " (" << a.end_gaps << " at the ends, "
       << (a.gaps - a.end_gaps) << " internal)\n"
       << "Identity: " << a.identity() << "% of all alignment columns\n";
  int overlap = a.length - a.end_gaps;
  if (overlap > 0)
    cout << "Overlap identity: " << 100.0 * a.matches / overlap << "% of the "
         << overlap << " columns between the first and last paired letters\n";

  if (shuffles > 0) {
    // Shuffling y keeps its length and letter composition but destroys its
    // order, so these scores show what an unrelated sequence would get.
    mt19937 rng(seed);
    string shuffled = y;
    vector<double> scores;
    int at_least = 0;
    for (int k = 0; k < shuffles; k++) {
      shuffle(shuffled.begin(), shuffled.end(), rng);
      int sc = nw_score_only(x, shuffled, s);
      scores.push_back(sc);
      if (sc >= a.score) at_least++;
    }
    double mean = accumulate(scores.begin(), scores.end(), 0.0) / shuffles;
    double var = 0;
    for (double v : scores) var += (v - mean) * (v - mean);
    double sd = sqrt(var / (shuffles > 1 ? shuffles - 1 : 1));

    cout << "\nShuffle test (" << shuffles << " shuffles of y, seed " << seed << "):\n"
         << "  Shuffled scores: mean " << mean << ", SD " << sd << "\n";
    if (sd > 0) cout << "  Real score is " << (a.score - mean) / sd << " SDs above the shuffled mean\n";
    cout << "  Shuffles scoring >= real score: " << at_least << " of " << shuffles << "\n";
  }
  return 0;
}
