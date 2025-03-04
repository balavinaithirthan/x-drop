//=== ydrop_one_sided_align ===

static score ydrop_one_sided_align
   (alignio*	io,
	int			reversed,			// true => search to the lower-left
	u8*			A,					// vertical sequence text
	u8*			B,					// horizontal sequence text
	unspos		M,					// vertical sequence length
	unspos		N,					// horizontal sequence length
	int			trimToPeak,
	editscript** script,
	unspos*		_end1,
	unspos*		_end2)
	{
	tback*		tb;					// tb is space provided to record traceback;
	u8*			tbp;				// .. tbp is the current traceback cell
	int			tbLen;				// .. and steps linearly from tbp->space;
	int			tbNeeded;			// .. tbRow[r] is the conceptual start of
									// .. the traceback cells for row r;  it is
									// .. indexed as tbRow[r][c] for c=LY..R',
									// .. where R' is the last cell used in
									// .. row r
	unspos		anchor1, anchor2;	// anchor positions in A and B
	scorerow*	allSub;				// substitution scores matrix
	score*		sub;				// substitution scores vector, current row
	score		gapOE, gapE;		// gap penalties
	score		yDrop;				// score drop threshold
	int			yDropTail;			// length of shortest score fall >= yDrop

									// nota bene:  row, col, leftCol, L, R,
									// .. LY, RY, prevLY, NN, npCol, end1
									// .. and end2 are all relative to the DP
									// .. matrix

	unspos		row, col = 0;		// current DP grid point
	unspos		leftCol;            // (copy of left column, for stats)
	sgnpos		L,  R;				// external column limits for current row
	unspos		LY, RY;				// actual column limits for current row;
									// .. ("Y" is for y-drop, not cartesian Y)
	unspos		prevLY;				// left column limit for previous row
	sgnpos		NN;					// truncated right side bound, current row
	unspos		npCol;				// last non-pruned cell in current row
	unspos		end1, end2;			// end of optimal alignment
	int			endIsBoundary;		// true => report boundaryScore instead of
									//         .. bestScore
	score		bestScore;			// score of best alignment seen so far
	score		boundaryScore;		// score of best alignment seen so far that
									// .. ends at the end of either sequence
	dpMatrix*	dynProg;			// DP cells
	dpCell*		dp;					// scans previous row of dp matrix
	dpCell*		dq;					// scans current row of dp matrix
	u8*			b;					// scans horizontal sequence
	score		c, d, i;			// running scores for DP cells
	score		cOpen, cNext, cTemp;// scratch values for cell scores
	activeseg*	active;				// list of segments that intersect the
									// .. sweep row within the feasible region
	galign*		alignList;
	galign*		rightAlign, *leftAlign;
	aliseg*		rightSeg,   *leftSeg;
	u8			link = 0;			// traceback link
	u8			op, prevOp;			// edit operations
	// sanity check;  if either sequence is empty there's no alignment

	if ((N <= 0) || (M <= 0))
		{ *(_end1) = *(_end2) = 0;  return 0; }

	// extract scoring constants

	allSub = io->scoring->sub;
	gapE   = io->scoring->gapExtend;
	gapOE  = io->scoring->gapOpen + gapE;	// (any new gap gets both penalties)
	yDrop  = io->yDrop;

	tb    = io->tb;
	tbLen = tb->size;

	if (gapE != 0)
		yDropTail = (yDrop/gapE) + 6;
	else
		{
		// when gapE is zero, the above results would be infinite;  but we can
		// limit yDropTail to the distance from the length of the sequence;
		// this can increase the amount of memory needed;  the solution here is
		// not completely sufficient;  "truncating alignment" reports are still
		// likely

		int maxYDropTail = 500*1000;
		if (N < (unsigned int) maxYDropTail) yDropTail = N+1;
		                                else yDropTail = maxYDropTail;
		}

	L = 0;
	R = N+1;								// (in blastz this was R=N)
	anchor1 = io->anchor1;
	anchor2 = io->anchor2;

	leftSeg = io->leftSeg;
	if (leftSeg != NULL)
		{
		L = signed_difference (leftSeg->b2, anchor2);
		if (leftSeg->type == diagSeg)
			L -= signed_difference (leftSeg->b1, anchor1);
		debugSnoopBounds_1;
		}

	rightSeg = io->rightSeg;
	if (rightSeg != NULL)
		{
		R = signed_difference (rightSeg->b2, anchor2);
		if (rightSeg->type == diagSeg)
			R -= signed_difference (rightSeg->b1, anchor1);
		debugSnoopBounds_2;
		}

	snoopSubprobsB_1;

	// if we're doing a reversed alignment we need to swap the L-R bounds (see
	// note (14))

	if (reversed)
		{
		sgnpos temp = 0;	// (placate compiler)
		if      ((leftSeg == NULL) && (rightSeg != NULL)) {               L = -R+1;  R = N+1;  }
		else if ((leftSeg != NULL) && (rightSeg == NULL)) { R    = -L-1;  L = 0;               }
		else if ((leftSeg != NULL) && (rightSeg != NULL)) { temp = -L-1;  L = -R+1;  R = temp; }
		debugSnoopBounds_3;
		}

	active     = NULL;
	rightAlign = io->rightAlign;
	leftAlign  = io->leftAlign;
	alignList  = (!reversed)? io->aboveList
	                        : io->belowList;

	// make sure we have a reasonable number of traceback rows to start with
	// (see note (3))

	tbrow_needed (minTbRowsNeeded);

	tbRow[0] = 0;
	tbp = tb->space;

	//////////
	// compute first row of dp matrix
	//////////

	// make sure we have enough traceback and dp space for the first row

	tbNeeded = yDropTail;
	if (tbNeeded > tbLen)
		suicide ("not enough space in trace_back array");

	dynProg = zalloc_or_die ("ydrop_one_sided_align dynProg", sizeof(dpMatrix));
	dp_ready (dynProg, tbNeeded);

	// compute first row of dp matrix (see notes (8), (10) and (13))

	dq = dynProg->p;
	dq->CC = cTemp = 0;						// set C[0][0]
	c = (dq++)->DD = -gapOE;				// set D[1][0]
	*(tbp++) = 0;

	snoopAlgorithm_2;

	for (col=1 ; (col<=N)&&(cTemp>=-yDrop) ; col++)
		{
		dq->CC = cTemp = c;					// set C[0][col]
		(dq++)->DD = c - gapOE;				// set D[1][col]
		c -= gapE;
		*(tbp++) = cFromI;

		snoopAlgorithm_3;
		}


	LY = 0;
	RY = col; // (1 column beyond the feasible region)


	end1 = end2 = 0;
	bestScore = 0;
	boundaryScore = negInf;
	endIsBoundary = false;

	for (row=1; row<=M ; row++)
		{

		snoopAlgorithm_4A;

		// update sweep row bounds, active segments, masking

		prevLY = LY;
		update_LR_bounds   (reversed,
		                    &rightSeg, &leftSeg, &rightAlign, &leftAlign,
		                    row, anchor1, anchor2, &L, &R, &LY, &RY);
		update_active_segs (reversed, &active, &alignList, dynProg->p-prevLY,
		                    row, anchor1, anchor2, LY, RY);

		snoopAlgorithm_4B;
		snoopSubprobsB_2;

		// make sure we have enough traceback and dp space for this row (see
		// note (3))

		tbrow_needed (row+1);
		gapped_extend_max_stat (maxDpRows, row+1);

		if (RY < LY) RY = LY;	// (see note 11)
		tbNeeded = RY - LY + yDropTail;
		if ((tbp - tb->space) + tbNeeded >= tbLen)
			{
			if (gapped_extend_inhibitTruncationReport)
				goto dp_finished;

			if (!reversed)
				fprintf (stderr, "truncating alignment ending at (" unsposCommaFmt ");",
				                 end1 + anchor1 + 1, end2 + anchor2 + 1);
			else
				fprintf (stderr, "truncating alignment starting at (" unsposCommaFmt ");",
				                 anchor1 + 2 - end1, anchor2 + 2 - end2);
			fprintf(stderr, "  anchor at (" unsposCommaFmt ")\n", anchor1, anchor2);

			if (!haveReportedTruncation)
				{
				haveReportedTruncation = true;
				fprintf(stderr, "truncation can be reduced by using --allocate:traceback to increase traceback memory\n");
				}
			goto dp_finished;
			}
		tbRow[row] = (tbp - tb->space) - LY;

		// set up DP pointers for this sweep row (see note (5))

		dp_ready (dynProg, tbNeeded);	// make sure we have enough DP space
		dq = dynProg->p;				// dq cells start at col == LY
		dp = dq + LY - prevLY;			// dp cells start at col == prevLY

		snoopAlgorithm_4C;

		sub = allSub[A[row]];

		col = leftCol = LY;
		b = B + col + 1;	// (b scans horizontal sequence, one column ahead)
		npCol = col;		// npCol records the last non-pruned position

		i = negInf;			// 'set' I[row][col]
		c = negInf;			// propose C[row][col]

		for ( ; (col<RY)&&((unspos)(b-B)<=N+1) ; col++)
			{

			d = dp->DD;						// get D[row][col]

			// at this point d, i and c contain the DP values for the cell at
			// [row][col];  the latter is the value for reaching C from previous
			// grid points, but since C also has edges entering it from the
			// current D and I nodes, we might improve it
			// nota bene: when we *can* improve c, we make an arbitrary choice
			// to prefer deletion to insertion (when i and d are equal)
			// nota bene 2: all paths through this series of ifs assign a value
			// to link

			if ((active != NULL) && (dp->mask == row))
				{ snoopAlgorithm_5;  prune;  snoopAlgorithm_5B;  continue; }

			if ((d > c) || (i > c))			// === we CAN improve C ===
				{
				// nota bene: both iExtend and dExtend are set here because
				// the value of the C and I (or C and D) are equal, so traceback
				// may as well take a gap extend into this cell
				if (d >= i) { c = d;  link = cFromD | iExtend | dExtend; }
				       else { c = i;  link = cFromI | iExtend | dExtend; }
				snoopAlgorithm_5;
				if (c < bestScore - yDrop)
					{ prune;  snoopAlgorithm_5B;  continue; }

#ifndef allowBackToBackGaps
				// not allowing back-to-back gaps, so we don't need to consider
				// opening a gap here
				i -= gapE;					// 'set' I[row][col+1]
				dq->DD = d - gapE;			// set D[row+1][col]
#else
				// back-to-back gaps are allowed, so we must consider gap opens
				cOpen = c - gapOE;
				d -= gapE;					// set D[row+1][col]
				if (cOpen > d) { dq->DD = cOpen;  link &= ~dExtend; }
				          else   dq->DD = d;

				i -= gapE;					// 'set' I[row][col+1]
				if (cOpen > i) { i = cOpen;  link &= ~iExtend; }
#endif // allowBackToBackGaps
				}
			else							// === we CANNOT improve C ===
				{
				snoopAlgorithm_5;
				if (c < bestScore - yDrop)
					{ prune;  snoopAlgorithm_5B;  continue; }

				if (c >= bestScore)
					{
					bestScore = c;  end1 = row;  end2 = col;  endIsBoundary = false;
					snoopAlgorithm_5A;
					}
				if ((!trimToPeak)
				      && (c >= boundaryScore)
				      && ((row == M) || (col == N)))
					{ boundaryScore = c;  end1 = row;  end2 = col;  endIsBoundary = true; }

				cOpen = c - gapOE;
				d -= gapE;					// set D[row+1][col]
				if (cOpen > d) { dq->DD = cOpen;  link = cFromC; }
				          else { dq->DD = d;      link = cFromC | dExtend; }

				i -= gapE;					// 'set' I[row][col+1]
				if (cOpen > i) i = cOpen;
				          else link |= iExtend;
				}

			npCol = col;					// save as last non-pruned position

			// save C for this column, and compute proposed C for the next
			// column (see note (6))

			cNext = (dp++)->CC+sub[*(b++)];	// propose C[row][col+1]
			(dq++)->CC = c;					// set C[row][col]
			c = cNext;
			*(tbp++) = link;				// set link into C[row][col]

			snoopAlgorithm_6;
			//snoopTraceback_1;
			}

		gapped_extend_add_stat (dpCellsVisited, col-leftCol);

		// if the feasible region is empty, we're done

		if (LY >= RY)
			goto dp_finished;

		// finish up this row, by either moving the right bound left or
		// prolonging the row to support an overhang on the row above

		snoopAlgorithm_7A;
		NN = ((rightSeg != NULL) && (R > 0))? (R-1) : ((sgnpos) N);

		if (RY > npCol+1)					// we hit ydrop prior to RY
			{
			RY = npCol+1;
			snoopAlgorithm_7B;
			}
		else
			{
			// the current row reached its right bound, but the row above may
			// still have a feasible overhang so we prolong this row with
			// insertions (see note (9))

			snoopAlgorithm_7C;
			while ((i >= bestScore - yDrop) && (((sgnpos)RY) <= NN))
				{
				if (((u32)(dq - dynProg->p)) >= dynProg->len)
					suicidef("(in ydrop_one_sided_align:%d, dq-dynProg->p==%d, dynProg->len=" unsposFmt ")",
					         __LINE__, dq - dynProg->p, dynProg->len);
				dq->CC = i;					// set C[row][col]
				(dq++)->DD = i - gapOE;		// set D[row+1][col]
				i -= gapE;					// 'set' I[row][col+1]
				*(tbp++) = cFromI;
				RY++;
				}
			snoopAlgorithm_7D;
			}

		// terminate the cell at the right boundary, so that nothing will
		// step from it (termination occurs if the if is false and we thus
		// *fail* to increment RY)

		if (((sgnpos)RY) <= NN)
			{
			if (((u32)(dq - dynProg->p)) >= dynProg->len)
				suicidef("(in ydrop_one_sided_align:%d, dq-dynProg->p==%d, dynProg->len=" unsposFmt ")",
				         __LINE__, dq - dynProg->p, dynProg->len);
			dq->DD = dq->CC = negInf;		// set D[row+1][col]
			RY++;							// .. and set C[row][col]
			snoopAlgorithm_7E;
			}
		}

dp_finished:
	snoopSubprobsB_3;
	*(_end1) = row = end1;
	*(_end2) = col = end2;

	snoopAlgorithm_8;

	//////////
	// traceback the alignment to create the edit script
	//////////

#ifdef snoopAlgorithm
	if (snoop)
		cTemp = 0;		// (place to set a breakpoint)
#endif // snoopAlgorithm

	for (prevOp=0 ; (row>=1) || (col>0) ; prevOp=op)
		{
		link = tb->space[tbRow[row] + col];
		op = link & cidBits;
		if ((prevOp == cFromI) && ((link & iExtend) != 0)) op = cFromI;
		if ((prevOp == cFromD) && ((link & dExtend) != 0)) op = cFromD;
		snoopTraceback_2

		if      (op == cFromI) {         col--;  edit_script_ins (script, 1); }
		else if (op == cFromD) { row--;          edit_script_del (script, 1); }
		else                   { row--;  col--;  edit_script_sub (script, 1); }
		snoopTraceback_3
		}

	filter_active_segs (&active, 2);	// (disposes of everything in the list)

	free_if_valid ("ydrop_one_sided_align dynProg->p", dynProg->p);
	free_if_valid ("ydrop_one_sided_align dynProg",    dynProg);

	if (endIsBoundary) return boundaryScore;
	              else return bestScore;
	}

//----------
//
// dp_ready--
//	Allocate the dynamic-programming "sweep row" array, and ensure that the
//	first n elements exist and have been initialized.
//
//----------
//
// Arguments:
//	dpMatrix*	dynProg:	The current sweep array.
//	unspos		needed:		The number of elements required.
//
// Returns:
//  nothing;  the contents of dynProg are (potentially) modified by adding
//	enough empty cells to satisfy the number needed
//
//----------
