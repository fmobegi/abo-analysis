#!/usr/bin/env python3


from Bio import pairwise2


class PairwiseAlignmentResult():
    
    
    def __init__(self, sequenceID, referenceSequence, querySequence):
        self.sequence1Aligned      = ''
        self.sequence2Aligned      = ''
        self.alignmentScore        = 0
        self.indexBegin            = 0
        self.indexEnd              = 0
        self.sequence1AlignedShort = ''
        self.sequence2AlignedShort = ''
        self.alignedSectionLength  = 0
        
        # Instead of storing a count of these events
        # I just want to store an array of the regions.
        # It's a tuple (begin, end)
        #self.matchTuples=[]
        #self.mismatchTuples=[]
        #self.insertionTuples=[]
        #self.deletionTuples=[]
        
        # Except for all-tuples. Lets just store these tuples:
        # (alignment state, begin, end)
        self.allTuples              = []
        
        #self.matchCount            = 0
        #self.mismatchCount         = 0
        #self.insertionCount        = 0
        #self.deletionCount         = 0
        #self.indelCount            = 0

        # SequenceIDs stores the ids of the alleles in this alignment.  
        self.sequenceID            = sequenceID
        self.sequence1             = referenceSequence
        self.sequence2             = querySequence
        
        self.alignPairwise()


    # TODO: Pairwise alignment takes a while.  
    # I can try some of the other alignment algorithms.
    # But then, why not just do a batch processing?


    def alignPairwise(self):      
        # Alignment Parameters.
        # A match score is 0.  I chose this so a perfect match is always a score of 0.
        # I chose to punish gaps more than mismatches.  
        # Likley there is tuning to be done here.  
        # I feel like match_score should be zero but I don't get any alignments when that happens.
        # Maybe the gap penalty should be the same as a mismatch penalty, because errors with repeats are somewhat common.

        match_score = 1
        mismatch_score = -1
        gap_open = -10
        gap_extend = -1


        # I want a local alignment, so it allows sequences of different lengths.  
        # The method to call depends on how i want to treat mismatches and indels.
        # I found some hints in the man pages for biopython's pairwise2 package
        # The match parameters are:

        #CODE  DESCRIPTION
        #x     No parameters. Identical characters have score of 1, otherwise 0.
        #m     A match score is the score of identical chars, otherwise mismatch
        #      score.
        #d     A dictionary returns the score of any pair of characters.
        #c     A callback function returns scores.
        #The gap penalty parameters are:

        #CODE  DESCRIPTION
        #x     No gap penalties.
        #s     Same open and extend gap penalties for both sequences.
        #d     The sequences have different open and extend gap penalties.
        #c     A callback function returns the gap penalties.

        # I chose localms.  So i can specify a mismatch score, and gap scores.
        # localms(sequenceA, sequenceB, match, mismatch, open, extend) -> alignments        
        #alignments = pairwise2.alignPairwise.localms(self.sequence1,self.sequence2
        #    ,match_score
        #    ,mismatch_score
        #    ,gap_open
        #    ,gap_extend)

        alignments = pairwise2.align.globalms(self.sequence1,self.sequence2
            ,match_score
            ,mismatch_score
            ,gap_open
            ,gap_extend)

        #Let's try a global alignment

        topAlignment = alignments[0]

        self.sequence1Aligned, self.sequence2Aligned, self.alignmentScore, self.indexBegin, self.indexEnd = topAlignment

        # I really don't need to store this information but why not?  Computer doesn't complain.
        self.sequence1AlignedShort = self.sequence1Aligned[self.indexBegin:self.indexEnd]
        self.sequence2AlignedShort = self.sequence2Aligned[self.indexBegin:self.indexEnd]

        #Count the number of mismatches etc.
        # An "insertion" is an extra base in sequence 2, not accounted in sequence 1.
        # A "deletion" is a base in sequence 1, missing in sequence 2
        # An indel is both of those.
        self.alignedSectionLength = len(self.sequence1AlignedShort)



        # Possible States = 'match' 'mismatch' 'insertion' 'deletion
        # Start at 0, because python enumerate method starts at 0.
        currentState = 'match'
        beginIndex = 0

        # Iterate the bases of the alignment, report regions that are match, mismatch, del, etc.
        for baseIndex, seq1Base in enumerate(self.sequence1AlignedShort):
            seq2Base = self.sequence2AlignedShort[baseIndex]
            if (seq1Base == seq2Base):
                if currentState != 'match':
                    self.storeAlignmentRegion(currentState, beginIndex, baseIndex-1)
                    beginIndex = baseIndex
                    currentState = 'match'

            elif (seq1Base == '-'):
                if currentState != 'insertion':
                    self.storeAlignmentRegion(currentState, beginIndex, baseIndex-1)
                    beginIndex = baseIndex
                    currentState = 'insertion'

            elif (seq2Base == '-'):
                if currentState != 'deletion':
                    self.storeAlignmentRegion(currentState, beginIndex, baseIndex-1)
                    beginIndex = baseIndex
                    currentState = 'deletion'

            elif currentState != 'mismatch':
                self.storeAlignmentRegion(currentState, beginIndex, baseIndex-1)
                beginIndex = baseIndex
                currentState = 'mismatch'

        # Guess I need to store the "last" aligned region.
        self.storeAlignmentRegion(currentState, beginIndex, baseIndex-1)
            
    def storeAlignmentRegion(self, state, beginIndex, endIndex):
        self.allTuples.append((state, beginIndex, endIndex))
        #if (state == 'match'):
        #    self.matchTuples.append((beginIndex, endIndex))
        #elif (state == 'mismatch'):
        #    self.mismatchTuples.append((beginIndex, endIndex))
        #elif (state == 'insertion'):
        #    self.insertionTuples.append((beginIndex, endIndex))
        #elif (state == 'deletion'):
        #    self.deletionTuples.append((beginIndex, endIndex))
        #else:
        #    raise Exception('I have a state that I do not know what to do with:' + str(state) + ',' + str(beginIndex) + ',' + str(endIndex))
            
