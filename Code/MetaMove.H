/*The goal of this class is to describe a graph edit move we make to generate a new graph.
 * We will make a move on a set of variables, typically a pair of variables, but may be more
 * if we are doing exact search. We need to store the type of move and the score
 * because of this new move. The move could be an edge addition, edge deletion, edge reversal,
 * Markov blanket addition, which is equal to addition of a set of edges */
#ifndef _MOVE_
#define _MOVE_

class MetaMove
{
	public:
		void setScoreImprovement(double);
		void setTargetMBScore(double);
		void setConditionSetInd(int);
		void setTFID(int);
		void setTargetID(int);
		int getTFID();
		int getTargetID();
		int getConditionSetInd();
		double getTargetMBScore();
		double getScoreImprovement();

	private:
		double scoreDelta;
		double targetMBScore;
		int regi;
		int targeti;
		int conditionSetInd;
};

#endif
