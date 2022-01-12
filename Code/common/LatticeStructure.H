/*The goal of this class is to store the Lattice structure of the random variables.
 * This class stores two types of relationships: superset relationship, subset relationship
 * Although, one can derive one from the other, we will store these relationships
 * in a matrix like data structure. For the superset relation, each row corresponds
 * to a factor and the entries in the row are factor ids that are contained in the factor
 * corresponding to this row. 
 * In the subset relation, each row corresponds to a factor and each entry in the row corresponds
 * to the factor ids that contain the factor corresponding to the row.
 * The superset relation is used to make the Markov blankets consistent, also while finding maximal 
 * clustes. 
 * The subset relation is used to search for Markov blankets
 * */
#ifndef _LATTICE_STRUCTURE
#define _LATTICE_STRUCTURE
#include "CommonTypes.H"

typedef map<int,INTINTMAP*> MATRIX;
typedef map<int,INTINTMAP*>::iterator MATRIX_ITER;

class LatticeStructure
{
	public:
		LatticeStructure();
		~LatticeStructure();
		int addSubset(int,int);
		int addSuperset(int,int);
		INTINTMAP* getSubsets(int);
		INTINTMAP* getSupersets(int);

		int getSupersets(int,INTINTMAP&,int);
		int getSuperSets(int fId,INTINTMAP& supersets,int currlevel,int level);

		int getAllSubsets(int, INTINTMAP&);

		int deleteFromLattice(int);
		int clear();
	private:
		MATRIX mSubset;
		MATRIX mSuperset;
};
#endif
