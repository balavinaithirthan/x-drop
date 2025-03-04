static score xdrop_extend_seed_hit
   (hitprocinfo*	hp,
	unspos*			_pos1,
	unspos*			_pos2,
	unspos*			_length)
	{
	unspos			pos1    = *_pos1;
	unspos			pos2    = *_pos2;
	unspos			length  = *_length;
	seq*			seq1    = hp->seq1;
	seq*			seq2    = hp->seq2;
	scoreset*		scoring = hp->scoring;
	score			xDrop   = hp->xDrop;
	sgnpos			diag, block2;
	unspos			oldDiagEnd, extent;
	u32				hDiag;
	u8*				s1, *s2, *stop, *leftStart, *rightStop, *rightBlock;
	score			similarity, runScore, leftScore, rightScore;
	int				adjustScore;

	//////////
	// extend to the left (loop 1)
	//
	// results:
	//	leftStart:  position of 1st nucleotide in extended match, in seq1
	//	leftScore:  score of bp added by left extension
	//	length:     possibly shortened (if the extension does not include the
	//	            .. entire hit)
	//////////

	s1 = seq1->v + pos1;	// start just past end of hit in both seq1 and seq2;
	s2 = seq2->v + pos2;	// .. we will pre-decrement before reads, so first
							// .. bp read are the ones at the tail of the hit

#ifdef extendHspFromLeft
	s1 = seq1->v + pos1 - length;
	s2 = seq2->v + pos2 - length;
#endif // extendHspFromLeft

	// determine stop location;  this is at the start of sequence 1, except
	// that if this diagonal ends (or is blocked) earlier in sequence 2, we
	// stop there
	// (see note 3;  instead of zero, we should use subsequence's start)

	if (unblockedLeftExtension) oldDiagEnd = 0;
	                       else oldDiagEnd = diagEnd[hDiag];
	block2 = (sgnpos) oldDiagEnd;
	if (block2 + diag > 0) stop = seq1->v + block2 + diag;
	                  else stop = seq1->v;

	// extend

	leftStart = s1;
	runScore  = leftScore = 0;

	while ((s1 > stop) && (runScore >= leftScore-xDrop))
		{
		snoopXDrop_Left;
		runScore += scoring->sub[*(--s1)][*(--s2)];
		if (runScore > leftScore)
			{
			leftStart = s1;
			leftScore = runScore;
			}
		}

	// adjust length if the extension is shorter than the hit

#ifndef extendHspFromLeft
	s2 = seq1->v + pos1 - length;	// (left end of hit)
	if (leftStart > s2)
		length -= leftStart - s2;
#endif // not extendHspFromLeft

#ifdef debugDiag
	p1 = s1;
	if (debugThisDiag)
		{
		p2 = leftStart;
		p3 = seq1->v+pos1-length;
		p4 = seq1->v+pos1;
		}
#endif
#ifdef collect_stats
	p1 = s1;
#endif

	//////////
	// extend to the right (loop 2)
	//
	// results:
	//	rightStop:  position of 1st nucleotide beyond extended match, in seq1
	//	similarity: increased by score of bp added by left and right extension
	//////////

	s1 = seq1->v + pos1;	// start just past end of hit in both seq1
	s2 = seq2->v + pos2;	// .. and seq2

#ifdef extendHspFromLeft
	s1 = seq1->v + pos1 - length;
	s2 = seq2->v + pos2 - length;
#endif // extendHspFromLeft

	// determine stop location;  this is at the end of sequence 1, except
	// that if this diagonal ends earlier in sequence 2, we stop there
	// (see note 3;  instead of sequence's end, we should use subsequence's end)

	block2 = (sgnpos) seq2->len;
	if ((sgnpos) seq1->len <= block2 + diag) stop = seq1->v + seq1->len;
	                                    else stop = seq1->v + block2 + diag;

	// extend

	rightStop = s1;
	runScore = rightScore = 0;

	while ((s1 < stop) && (runScore >= rightScore-xDrop))
		{
		snoopXDrop_Right;
		runScore += scoring->sub[*(s1++)][*(s2++)];
		if (runScore > rightScore)
			{
			rightStop  = s1;
			rightScore = runScore;
			}
		}
	rightBlock = s1;

	// adjust length if the extension is shorter than the hit

#ifdef extendHspFromLeft
	s2 = seq1->v + pos1;			// (past right end of hit)
	if (rightStop < s2)
		length -= s2 - rightStop;
#endif // extendHspFromLeft

#ifdef debugDiag
	if (debugThisDiag)
		dump_extended_match (stdout, seq1, seq2, diag,
		                     p1, p2, p3, p4, rightStop, s1);
#endif

	similarity = leftScore + rightScore;

	//////////
	// (for debugging only)
	// determine if some subrange of the HSP outscores the whole HSP
	//
	// We use the algorithm described in Bentley's "Programming Pearls" ("A
	// scanning Algorithm", section 8.4, page 81 in the second edition).
	//
	// bestSubrangeScore == Bentley's maxSoFar
	// subrangeScore     == Bentley's maxEndingHere
	//////////

#ifdef snoopHspSubrange

	{
	score subrangeScore;
	u8*   currLeft, *subLeft, *subRight;

	s1 = leftStart;
	s2 = seq2->v + diagToPos2 (diag, leftStart - seq1->v);
	currLeft = s1;
	subLeft = subRight = s1;

	subrangeScore = bestSubrangeScore = 0;
	runScore = 0;
	while (s1 < rightStop)
		{
		runScore      = runScore      + scoring->sub[*s1][*s2];
		subrangeScore = subrangeScore + scoring->sub[*(s1++)][*(s2++)];
		if (subrangeScore < 0)
			{ subrangeScore = 0;  currLeft = s1; }
		if (subrangeScore > bestSubrangeScore)
			{
			bestSubrangeScore = subrangeScore;
			subLeft  = currLeft;
			subRight = s1;
			}

		if (debugThisDiag)
			printf (unsposFmt ":"
			        " %c%c " scoreFmtSimple
			        " " scoreFmtSimple
			        " " scoreFmtSimple " " unsposFmt
			        " " scoreFmtSimple " " unsposFmt ".." unsposFmt
			        "\n",
			        (s1-1) - seq1->v,
			        s1[-1], s2[-1], scoring->sub[s1[-1]][s2[-1]],
			        runScore,
			        subrangeScore, currLeft - seq1->v,
			        bestSubrangeScore, subLeft - seq1->v, subRight - seq1->v);
		} 

	hspPos1 = hspPos2 = hspLen = subPos1 = subPos2 = subLen = 0;

	if (bestSubrangeScore > similarity)
		{
		hspPos1 = leftStart - seq1->v;
		hspPos2 = diagToPos2 (diag, leftStart - seq1->v);
		hspLen  = rightStop - leftStart;
		subPos1 = subLeft - seq1->v;
		subPos2 = diagToPos2 (diag, subLeft - seq1->v);
		subLen  = subRight - subLeft;
		seed_search_count_stat (suboptimalHsp);
		}
	}

#endif // snoopHspSubrange

	//////////
	// record the extent of HSP search on this diagonal
	//////////

	// record the extent

	extent = (unspos) (((sgnpos) (rightBlock-seq1->v)) - diag);
	if (extent > diagEnd[hDiag])
		{
		diagEnd   [hDiag] = extent;
		diagActual[hDiag] = diag;
#ifdef snoopDiagHash
		fprintf (stderr, "  setting    diag %9s"
		                 ", diagActual[%04X] = %9s"
		                 ", diagEnd[%04X] = " unsposFmt
		                 ", seed end = " unsposSlashFmt
		                 ", in seq 1: " unsposDotsFmt "\n",
		                 pair_diagonal_as_text(pos1,pos2),
		                 hDiag, diagonal_as_text(diagActual[hDiag]),
		                 hDiag, diagEnd[hDiag],
		                 pos1, pos2, start2, pos2);
#endif // snoopDiagHash
		}

#ifdef debugDiag
	if (debugThisDiag)
		{
		printf ("gfex: (diag %9s)      " unsposSlashFmt " diagEnd[%04X] <-- " unsposFmt,
				pair_diagonal_as_text(pos1,pos2), pos1, pos2,
				hDiag, diagEnd[hDiag]);
		printf (" (" unsposSlashFmt " -> " unsposSlashFmt ")\n",
				(unspos) (rightStop-seq1->v),  (unspos) (rightStop-seq1->v-diag),
				(unspos) (rightBlock-seq1->v), diagEnd[hDiag]);
		}
#endif

	snoopXDrop_Score;

#ifdef collect_stats
	seed_search_add_stat (bpExtended, rightBlock-p1);
#endif

	//////////
	// update length of hit
	//////////

	pos1   = (unspos) (rightStop - seq1->v);
	pos2   = (unspos) (((sgnpos) pos1) - diag);
	length = (unspos) (rightStop - leftStart);

	//////////
	// if the extended hit's score is acceptable, but not very high, adjust
	// the score downward, based on the entropy of the sequences in the match;
	// note that we only adjust positive scores (since hspZeroThreshold ==
	// max(0,hspZeroThreshold)), otherwise the entropy adjustment would
	// *increase* the score when entropy is poor.
	//
	// When an adaptive scoring threshold is being used, we can't determine
	// what a reasonable "high enough" threshold is to not bother to perform
	// the entropy reduction, so we perform the reduction on any extended hit
	// that could potentially make the hsp table.
	//
	// $$$ Heuristically, we could still estimate some reasonable "high enough"
	//     .. threshold based on the current low score threshold, how many hits
	//     .. we have accepted/rejected so far, and how much room is left in
	//     .. the table.  Simpler schemes may also work well in practice.  This
	//     .. should be considered if the entropy calculation ends up being a
	//     .. significant time factor.
	//////////

	// decide whether to adjust the score

	if (!hp->entropicHsp)
		adjustScore = false;
	else if (hp->hspThreshold.t == 'S')	// (fixed score threshold)
		adjustScore = (similarity >=   hp->hspZeroThreshold)
	               && (similarity <= 3*hp->hspThreshold.s);
	else if (similarity <= 0)			// (adaptive score threshold, negative)
		adjustScore = false;
	else								// (adaptive score threshold, positive)
		{
		segtable* anchors = *(hp->anchors);
		adjustScore = (anchors->len > 0)
		           && (similarity >= anchors->lowScore);
		}

	// adjust it

	if (adjustScore)
		{
		double q = entropy (hp->seq1->v + pos1 - length,
		                    hp->seq2->v + pos2 - length,
		                    length);

		score rawS = similarity;
		similarity *= q;
		if ((similarity < hp->hspThreshold.s) && (hp->reportEntropy))
			fprintf(stderr, "hit of score " scoreFmtSimple
			                " at " unsposSlashFmt "#" unsposFmt " (diag " sgnposFmt " had block at " unsposFmt ")"
			                " fails entropy filter (%f)\n",
							rawS,
							pos1-length, pos2-length, length,
							diag, oldDiagEnd,
							q);
#ifdef snoopEntropy
		else
			fprintf(stderr, "hit of score " scoreFmtSimple
			                " at " unsposSlashFmt "#" unsposFmt " (diag " sgnposFmt " had block at " unsposFmt ")"
			                " passes entropy filter (%f)\n",
							rawS,
							pos1-length, pos2-length, length,
							diag, oldDiagEnd,
							q);
#endif // snoopEntropy

#ifdef snoopHspSubrange
		bestSubrangeScore *= q;
#endif // snoopHspSubrange
		}

	//////////
	// decide whether or not this extended seed hit is an hsp.
	//////////

	dbg_timing_count_stat (ungappedExtensions);

	// if it doesn't score high enough, discard it

	if ((hp->hspThreshold.t == 'S')				// (fixed score threshold)
	 && (similarity < hp->hspThreshold.s))
		{
		seed_search_count_stat (lowScoringHsps);
#ifdef snoopHspSubrange
		if (bestSubrangeScore >= hp->hspThreshold.s)
			{
			fprintf (stderr, "WARNING: discarded HSP " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple
							 " but subrange " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple "\n",
			                 hspPos1, hspPos2, hspLen, similarity,
			                 subPos1, subPos2, subLen, bestSubrangeScore);

			seed_search_count_stat (suboptimalHspB);
			}
		else if ((seed_search_dbgSubrangeHsps)
		      && (bestSubrangeScore > similarity))
			fprintf (stderr, "INFO: HSP " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple
							 " but subrange " unsposSlashFmt "#" unsposFmt
							 " scores " scoreFmtSimple "\n",
							 hspPos1, hspPos2, hspLen, similarity,
							 subPos1, subPos2, subLen, bestSubrangeScore);
#endif // snoopHspSubrange
		return noScore;
		}

#ifdef snoopHspSubrange
	if ((seed_search_dbgSubrangeHsps)
	 && (bestSubrangeScore > similarity))
		fprintf (stderr, "INFO: HSP " unsposSlashFmt "#" unsposFmt
						 " scores " scoreFmtSimple
						 " but subrange " unsposSlashFmt "#" unsposFmt
						 " scores " scoreFmtSimple "\n",
						 hspPos1, hspPos2, hspLen, similarity,
						 subPos1, subPos2, subLen, bestSubrangeScore);
#endif // snoopHspSubrange

	// it's a keeper

	*_pos1   = pos1;
	*_pos2   = pos2;
	*_length = length;

	if (hp->anchors != NULL)
		(*(hp->anchors))->haveScores = true;

	seed_search_count_stat (hsps);
	dbg_timing_count_stat  (hsps);

	return similarity;
	}