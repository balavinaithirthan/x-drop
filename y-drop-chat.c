/***********************************************************************
 * A BASIC ONE-SIDED Y-DROP ALIGNMENT EXAMPLE
 *
 * Compile & run (e.g.):
 *   gcc -o ydrop_example ydrop_example.c
 *   ./ydrop_example
 *
 * This code is simplified from LASTZ's ydrop_one_sided_align() -- all
 * the advanced bounds logic, pruning, debug printouts, and segment‐list
 * logic is omitted to keep things straightforward.
 ***********************************************************************/

 #include <stdio.h>
 #include <stdlib.h>
 #include <limits.h>
 
 /*----------------------------------------------------------------------*
  * Basic type and constant definitions
  *----------------------------------------------------------------------*/
 
 typedef unsigned int  unspos;   // for sequence positions
 typedef int           score;    // for DP cell scores
 
 // A convenient representation of negative infinity:
 #define NEG_INF  (INT_MIN / 2)
 
 // Bits for traceback encoding:
 #define cFromC   0x00  // came from Diagonal (substitution)
 #define cFromD   0x01  // came from Deletion
 #define cFromI   0x02  // came from Insertion
 #define dExtend  0x04  // extended a deletion gap
 #define iExtend  0x08  // extended an insertion gap
 
 /***********************************************************************
  * Simple data structures to illustrate the concept
  ***********************************************************************/
 
 // A “scoring” structure with minimal info: gap penalties + substitution table
 typedef struct {
     score gapOpen;          // penalty for opening a gap
     score gapExtend;        // penalty per gap extension
     score sub[256][256];    // a small substitution matrix (indexed by residue)
 } scoring_t;
 
 // A structure that carries relevant data for the alignment call
 typedef struct {
     scoring_t* scoring;     // pointer to scoring parameters
     score      yDrop;       // y‐drop threshold
 } alignio;
 
 // Each DP cell needs to keep track of at least the best "match/mismatch" path
 // (C) and the best "deletion" path (D); we will keep “insertion” path (I)
 // in a local variable during the row sweep.
 typedef struct {
     score CC;  // best alignment score ending with a match or mismatch
     score DD;  // best alignment score ending with a deletion gap
 } dpCell;
 
 // Our DP “matrix” for the y-drop row-sweep approach. We only store one row at
 // a time for the previous row and one row for the current row (but for clarity
 // we put them together in a single array here).
 typedef struct {
     dpCell* p;     // the DP cells
     size_t  len;   // how many dpCell's are currently allocated
 } dpMatrix;
 
 // Simple container for traceback bits
 typedef struct {
     unsigned char* space;
     size_t         size;
 } tback;
 
 // Minimal edit-script type, and stubs for writing to it
 typedef struct {
     // for demonstration, you could store the actual ops here
     // but we leave it empty
 } editscript;
 
 static void edit_script_ins(editscript** script, int length) { /* stub */ }
 static void edit_script_del(editscript** script, int length) { /* stub */ }
 static void edit_script_sub(editscript** script, int length) { /* stub */ }
 
 /***********************************************************************
  * Utility function to ensure we have enough DP cells allocated
  ***********************************************************************/
 static void dp_ready(dpMatrix* dp, size_t needed)
 {
     // If we don't have enough cells in dp->p, expand it
     if (dp->len < needed) {
         dp->p = (dpCell*)realloc(dp->p, sizeof(dpCell)*needed);
         if (dp->p == NULL) {
             fprintf(stderr, "Out of memory in dp_ready().\n");
             exit(EXIT_FAILURE);
         }
         dp->len = needed;
     }
 }
 
 /***********************************************************************
  * The core y-drop DP function (a simplified version)
  *
  *   io         -- carries gap penalties, sub matrix, yDrop
  *   reversed   -- if non-zero, we'd do the DP 'backwards' (omitted here)
  *   A, B       -- sequences to align (assume zero-based in this example)
  *   M, N       -- lengths of A and B
  *   trimToPeak -- if non-zero, might trim alignment at local maximum (omitted)
  *   script     -- pointer to the alignment edit script (for traceback)
  *   _end1,_end2-- where the alignment ends in A,B (outputs)
  *
  * Return value -- the best alignment score found
  ***********************************************************************/
 score ydrop_one_sided_align
 (
     alignio*    io,
     int         reversed,   // (ignored in this basic example)
     unsigned char* A,       // vertical sequence
     unsigned char* B,       // horizontal sequence
     unspos      M,          // length of A
     unspos      N,          // length of B
     int         trimToPeak, // (ignored here)
     editscript** script,
     unspos*     _end1,
     unspos*     _end2
 )
 {
     // if either sequence is empty, no alignment
     if ((M == 0) || (N == 0)) {
         *_end1 = 0;
         *_end2 = 0;
         return 0;
     }
 
     // gather scoring constants
     score gapE   = io->scoring->gapExtend;
     score gapOE  = io->scoring->gapOpen + gapE; // cost when opening a gap
     score yDrop  = io->yDrop;
 
     // allocate a big block for traceback bits
     // In a full y-drop approach, we’d need something more dynamic;
     // we do a naive approach for demonstration
     tback tb;
     tb.size = (M+1) * (N+1);           // large enough for this example
     tb.space = (unsigned char*)malloc(tb.size);
     if (tb.space == NULL) {
         fprintf(stderr, "Out of memory for traceback.\n");
         exit(EXIT_FAILURE);
     }
     // We'll fill tb.space row-by-row with our DP transitions
     //   (tb.space[row*(N+1) + col] = link)
     // link is a combination of cFromX plus any dExtend/iExtend bits
 
     // allocate DP matrix
     dpMatrix dynProg;
     dynProg.p = NULL;
     dynProg.len = 0;
 
     // Make sure we have enough dpCell's for at least 2 rows of length (N+1),
     // plus some margin. We'll just allocate (M+1)*(N+1) to keep it simple.
     dp_ready(&dynProg, (M+1)*(N+1));
     // We’ll treat row i as &dynProg.p[ i*(N+1) ].
 
     // Initialize the first row (row=0)
     //   In standard alignment, you'd have a typical DP boundary condition,
     //   e.g. scoring a gap from col=0..N
     {
         dpCell* row0 = &dynProg.p[0*(N+1)];
         row0[0].CC = 0;       // match portion at (0,0) is 0
         row0[0].DD = -gapOE;  // “deletion” at next row’s row=1, col=0
         tb.space[0*(N+1) + 0] = 0;  // link bits for (0,0)
 
         // C(0,col) extends insertion from left
         // For demonstration, we’ll do a typical linear “infinite gap penalty”
         score c = row0[0].DD; // start from D(1,0) idea
         for (unspos col=1; col <= N; col++) {
             row0[col].CC = c;           // continuing insertion
             row0[col].DD = c - gapE;    // next row’s deletion
             tb.space[0*(N+1) + col] = cFromI;
             c -= gapE;
         }
     }
 
     /*******************************************************************
      * The main y-drop row-sweep
      *******************************************************************/
     score bestScore = 0;
     unspos bestRow = 0, bestCol = 0;
 
     for (unspos row = 1; row <= M; row++)
     {
         // Pointers to the previous row and current row in DP storage
         dpCell* prev = &dynProg.p[(row-1)*(N+1)];
         dpCell* curr = &dynProg.p[ row   *(N+1)];
 
         // We'll track the insertion path in a local variable
         // (like I[row][col]) while we sweep from col=0..N
         score i = NEG_INF;
 
         // Grab the substitution row for the residue A[row-1]
         // (assuming A and B are 0-based arrays)
         unsigned char aRes = A[row-1];
         const score* subRow = io->scoring->sub[aRes];
 
         // At col=0, we typically have:
         curr[0].CC = NEG_INF;                  // can't align anything at col=0
         curr[0].DD = prev[0].CC - gapOE;       // open a deletion from the top
         tb.space[row*(N+1) + 0] = cFromD;      // marks that we came from above
 
         // Walk across columns
         for (unspos col = 1; col <= N; col++)
         {
             // D = best path that leads to a deletion at (row,col)
             score d = prev[col].DD - gapE;  // extends a deletion from above
             // or open new deletion from curr[col].CC
             score cOpen = (curr[col].CC = NEG_INF); // we will fill CC soon
             // We'll fill cOpen after we get CC for this cell
 
             // First compute the diagonal “match/mismatch” from (row-1,col-1).
             // That is prev[col-1].CC plus sub cost from A[row-1], B[col-1].
             score cVal = prev[col-1].CC + subRow[B[col-1]];
 
             // We see if either D or I can improve cVal
             if (d > cVal || i > cVal)
             {
                 if (d >= i) {
                     cVal = d;
                     tb.space[row*(N+1) + col] = cFromD | dExtend | iExtend;
                 } else {
                     cVal = i;
                     tb.space[row*(N+1) + col] = cFromI | dExtend | iExtend;
                 }
             }
             else
             {
                 // we came from the diagonal
                 tb.space[row*(N+1) + col] = cFromC;
             }
 
             // The final CC is cVal
             curr[col].CC = cVal;
 
             // Update best score if needed
             if (cVal > bestScore) {
                 bestScore = cVal;
                 bestRow   = row;
                 bestCol   = col;
             }
 
             // Once we know cVal == C[row][col], we can define cOpen
             // to open a gap from this position
             cOpen = cVal - gapOE;
 
             // Now set D[row+1][col], i.e. “DD” for the current row
             {
                 score fromD = prev[col].DD - gapE; // continuing a deletion
                 if (cOpen > fromD) {
                     curr[col].DD = cOpen; // open deletion from cVal
                     // remove dExtend bit if you like
                 } else {
                     curr[col].DD = fromD; // extend deletion from above
                 }
             }
 
             // Finally update I[row][col+1] by either opening from cVal or
             // extending from i
             {
                 score fromI = i - gapE;
                 if (cOpen > fromI) {
                     i = cOpen; // open insertion from cVal
                 } else {
                     i = fromI; // extend insertion from left
                 }
             }
 
             // Y-drop check: if cVal < bestScore - yDrop,
             // you might prune here in a more complete implementation.
         }
     }
 
     // The best cell is (bestRow, bestCol) with bestScore
     *_end1 = bestRow;
     *_end2 = bestCol;
 
     /*******************************************************************
      * Traceback to build the edit script
      *******************************************************************/
     {
         unspos row = bestRow;
         unspos col = bestCol;
         unsigned char prevOp = cFromC;  // arbitrary start
 
         while (row > 0 || col > 0)
         {
             unsigned char link = tb.space[row*(N+1) + col];
             unsigned char op   = link & 0x03; // cFromI, cFromD, or cFromC
 
             // If previous op was I and iExtend is set, continue that insertion
             if ((prevOp == cFromI) && (link & iExtend)) op = cFromI;
             // If previous op was D and dExtend is set, continue that deletion
             if ((prevOp == cFromD) && (link & dExtend)) op = cFromD;
 
             // Go whichever direction the op indicates
             if (op == cFromI) {
                 // insertion => move left
                 col--;
                 edit_script_ins(script, 1);
             }
             else if (op == cFromD) {
                 // deletion => move up
                 row--;
                 edit_script_del(script, 1);
             }
             else {
                 // diagonal => match/mismatch => move diagonally
                 row--;
                 col--;
                 edit_script_sub(script, 1);
             }
             prevOp = op;
         }
     }
 
     // Clean up
     free(tb.space);
     free(dynProg.p);
 
     return bestScore;
 }
 
 /***********************************************************************
  * A SIMPLE TEST HARNESS
  *----------------------------------------------------------------------*/
 int main(void)
 {
     // Example sequences
     // For clarity, we pretend 1..n are valid residues
     unsigned char A[] = {1,2,3,2,1};
     unsigned char B[] = {1,2,2,3,1,2};
     unspos M = 5;
     unspos N = 6;
 
     // Create a simple scoring table: +1 for match, -1 for mismatch
     scoring_t scoring;
     for (int x = 0; x < 256; x++) {
         for (int y = 0; y < 256; y++) {
             scoring.sub[x][y] = (x == y)? 1 : -1;
         }
     }
     scoring.gapOpen   = 2;
     scoring.gapExtend = 1;
 
     // Make alignio with a y-drop threshold
     alignio io;
     io.scoring = &scoring;
     io.yDrop   = 5;   // e.g. drop out once score is 5 below best
 
     // Prepare an edit script pointer (we don’t really store ops here)
     editscript* script = NULL;
 
     // End positions
     unspos end1, end2;
 
     // Perform the alignment
     score best = ydrop_one_sided_align(
         &io,
         0,           // reversed == false
         A,
         B,
         M,
         N,
         0,           // trimToPeak == false
         &script,
         &end1,
         &end2
     );
 
     // Report
     printf("Best score = %d\n", best);
     printf("Alignment ends at A:%u, B:%u\n", end1, end2);
 
     return 0;
 }
 