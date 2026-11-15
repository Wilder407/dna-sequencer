// Needleman-Wunsch global sequence alignment
// Author: Lisa Wilder (GitHub: Wilder407)

#ifndef MATRIX_H__
#define MATRIX_H__

#include <string>
#include <vector>

// Scores used to fill the matrix.
//   match:    added when two aligned letters are the same
//   mismatch: added when two aligned letters differ
//   gap:      added when a letter is aligned to a gap (an insertion or deletion)
struct Scoring {
  int match = 1;
  int mismatch = -1;
  int gap = -2;
  // Semi-global mode: gaps at the start or end of either sequence are free.
  // Use when one sequence is a partial read (for example an EST) of the other.
  bool free_end_gaps = false;
};

// Direction the best score came from, used for the traceback.
enum class Dir : char { None = ' ', Diag = '\\', Up = '|', Left = '-' };

struct NWMatrices {
  std::vector<std::vector<int>> score;  // F, size (m+1) x (n+1)
  std::vector<std::vector<Dir>> dir;    // traceback directions, same size
};

struct Alignment {
  int score = 0;           // optimal global alignment score, F[m][n]
  std::string aligned_x;   // sequence x with '-' for gaps
  std::string aligned_y;   // sequence y with '-' for gaps
  int matches = 0;
  int mismatches = 0;
  int gaps = 0;            // number of gap positions (internal + end)
  int end_gaps = 0;        // gap positions at either end of the alignment
  int length = 0;          // alignment length (columns)
  double identity() const { return length ? 100.0 * matches / length : 0.0; }
};

// Cleans a raw sequence: drops FASTA/comment header lines (starting with '>'
// or '//'), whitespace and digits, and converts to uppercase.
std::string clean_sequence(const std::string& raw);

// Reads every sequence in a FASTA-style file. Records start with a line
// beginning with '>' or '//'. A file with no headers is one sequence.
std::vector<std::string> read_sequences(const std::string& path);

// True if seq is non-empty and contains only A, C, G, T, U (uppercase).
bool valid_input(const std::string& seq);

int match_score(char a, char b, const Scoring& s);

// Fills the score matrix F and the traceback direction matrix.
NWMatrices N_W(const std::string& x, const std::string& y, const Scoring& s);

// Score-only version using two rows of memory, for long sequences and
// permutation tests where the traceback is not needed.
int nw_score_only(const std::string& x, const std::string& y, const Scoring& s);

// Follows the traceback to F[0][0] to build the alignment. The alignment
// score is F[m][n] (global) or the best cell in the last row or column
// (free end gaps); the traceback recovers which letters pair.
Alignment traceback(const NWMatrices& M, const std::string& x, const std::string& y,
                    const Scoring& s);

// Prints the score and direction matrices with sequence labels.
void print_matrix(const NWMatrices& M, const std::string& x, const std::string& y);

// Prints the alignment in blocks of `width` columns with a match line
// ('|' match, '.' mismatch, ' ' gap).
void print_alignment(const Alignment& a, int width = 60);

#endif  // MATRIX_H__
