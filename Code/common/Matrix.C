
#include"Matrix.H"
#include <math.h>
#include <map>
#include <vector>
#include <chrono>
using namespace std::chrono;


Matrix::Matrix()
{
	matrix=NULL;
	row=0;
	col=0;
}

Matrix::Matrix(int r,int c) : 
	row(r),
	col(c)
{
	matrix=gsl_matrix_alloc(row,col);
}

Matrix::~Matrix()
{
	if(matrix!=NULL)
	{
		gsl_matrix_free(matrix);
		matrix=NULL;
	}
}

int 
Matrix::init(int r,int c)
{
	row=r;
	col=c;
	matrix=gsl_matrix_alloc(r,c);
	return 0;
}

int
Matrix::initAsIdentity()
{
	if((row==0)||(col==0))
	{
		cout << "Warning!! Row or col = 0" << endl;
	}
	gsl_matrix_set_identity (matrix);
	return 0;
}
	
//Add matrix b to this matrix
//Return NULL on failure
//Caller must free this matrix
Matrix* 
Matrix::addMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "Dimensions do not match" << endl;
		return NULL;	
	}
	Matrix* res=new Matrix(row,col);
	gsl_matrix_memcpy (res->matrix, matrix);
	gsl_matrix_add(res->matrix,b->matrix);
	return res;	
}
Matrix* 
Matrix::subtractMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "Dimensions do not match" << endl;
		return NULL;	
	}
	Matrix* res=new Matrix(row,col);
	gsl_matrix_memcpy (res->matrix, matrix);
	gsl_matrix_sub(res->matrix,b->matrix);
	return res;	
}
/*
//Multiply this with b
Matrix* 
Matrix::multiplyMatrix(Matrix* b)
{
	if(col!=b->getRowCnt())
	{
		return NULL;
	}
	Matrix* res=new Matrix(row,b->getColCnt());
	gsl_matrix_set_zero (res->matrix);

	gsl_matrix_float* resmatrix=gsl_matrix_float_alloc(row,b->getColCnt());
	convertToFloat(resmatrix,res->matrix,row,b->getColCnt());

	gsl_matrix_float* cmatrix=gsl_matrix_float_alloc(row,col);
	convertToFloat(cmatrix,matrix,row,col);

	gsl_matrix_float* bmatrix=gsl_matrix_float_alloc(b->getRowCnt(),b->getColCnt());
	convertToFloat(bmatrix,b->matrix,b->getRowCnt(),b->getColCnt());
	
	gsl_blas_sgemm (CblasNoTrans, CblasNoTrans, 1, cmatrix, bmatrix, 0, resmatrix);
	convertFromFloat(resmatrix,res->matrix,row,b->getColCnt());
	gsl_matrix_float_free(resmatrix);
	gsl_matrix_float_free(cmatrix);
	gsl_matrix_float_free(bmatrix);
	return res;	
}*/

Matrix*
Matrix::multiplyMatrix(Matrix* b)
{
    if(col!=b->getRowCnt())
    {
        return NULL;
    }
    Matrix* res=new Matrix(row,b->getColCnt());
    gsl_matrix_set_zero(res->matrix);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, matrix, b->matrix, 0, res->matrix);
    return res;
}


int 
Matrix::addWithMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "Dimensions do not match" << endl;
		return -1;	
	}
	gsl_matrix_add(matrix,b->matrix);
	return 0;	
}
	
int 
Matrix::subtractWithMatrix(Matrix* b)
{
	if(!dimequal(b))
	{
		cout << "Dimensions do not match" << endl;
		return -1;	
	}
	gsl_matrix_sub(matrix,b->matrix);
	return 0;	
}

// no need to transform to float, matrix is modified
int
Matrix::multiplyWithMatrix(Matrix* b)
{
    if(col!=b->getRowCnt())
    {
        return -1;
    }
    gsl_matrix* res=gsl_matrix_alloc(row,b->getColCnt());
    gsl_matrix_set_zero(res);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, matrix, b->matrix, 0, res);
    gsl_matrix_memcpy(matrix, res);
    gsl_matrix_free(res);
    return 0;
}
/*
int 
Matrix::multiplyWithMatrix(Matrix* b)
{
	if(col!=b->getRowCnt())
	{
		return -1;
	}
	
	gsl_matrix* res=gsl_matrix_alloc(row,b->getColCnt());
	gsl_matrix_set_zero(res);
	
	gsl_matrix_float* resmatrix=gsl_matrix_float_alloc(row,b->getColCnt());
	convertToFloat(resmatrix,res,row,b->getColCnt());

	gsl_matrix_float* cmatrix=gsl_matrix_float_alloc(row,col);
	convertToFloat(cmatrix,matrix,row,col);

	gsl_matrix_float* bmatrix=gsl_matrix_float_alloc(b->getRowCnt(),b->getColCnt());
	convertToFloat(bmatrix,b->matrix,b->getRowCnt(),b->getColCnt());
	
	gsl_blas_sgemm(CblasNoTrans, CblasNoTrans, 1, cmatrix, bmatrix, 0, resmatrix);
	convertFromFloat(resmatrix,res,row,b->getColCnt());

	gsl_matrix_float_free(resmatrix);
	gsl_matrix_float_free(cmatrix);
	gsl_matrix_float_free(bmatrix);
	
	gsl_matrix_memcpy(matrix, res);
	gsl_matrix_free(res);
	return 0;
}*/
	
int 
Matrix::addScalar(double aVal)
{
	gsl_matrix_add_constant(matrix,aVal);
	return 0;
}

int 
Matrix::subtractScalar(double aVal)
{
	gsl_matrix_add_constant(matrix,-1*aVal);
	return 0;
}

int 
Matrix::multiplyScalar(double aVal)
{
	gsl_matrix_scale(matrix,aVal);
	return 0;
}

int 
Matrix::divideScalar(double aVal)
{
	if(aVal==0)
	{
		return -1;
	}
	gsl_matrix_scale(matrix,1/aVal);
	return 0;
}

int 
Matrix::setValue(double val,int i,int j)
{
	gsl_matrix_set(matrix, i, j, val);
	return 0;
}

int 
Matrix::setAllValues(double val)
{
	gsl_matrix_set_all (matrix,val);
	return 0;
}

double
Matrix::getValue(int i,int j)
{
	double val=gsl_matrix_get(matrix, i, j);
	return val;
}

// extract submatrix using index:
Matrix*
Matrix::getSubMatrix(vector<int>& idx)
{
    int subrow=idx.size();
    gsl_matrix* smat=gsl_matrix_alloc(subrow,col);
    for(int i=0;i<subrow;i++)
    {
        gsl_vector * x1 = gsl_vector_alloc(col);
        gsl_matrix_get_row(x1, matrix, idx[i]);
        gsl_matrix_set_row(smat,i,x1);
        gsl_vector_free(x1);
    }
    Matrix* submat=new Matrix(subrow,subrow);
    //gsl_matrix* submat=gsl_matrix_alloc(subrow,subrow);
    // get columns
    for(int j=0;j<subrow;j++)
    {
        gsl_vector * x1 = gsl_vector_alloc(subrow);
        gsl_matrix_get_col(x1, smat, idx[j]);
        gsl_matrix_set_col(submat->matrix,j,x1);
        gsl_vector_free(x1);
    }
    gsl_matrix_free(smat);
    //submat->showMatrix(cout);
    return submat;
}

Matrix*
Matrix::getSubMatrix(int k1, int k2, int n1, int n2)  // get from 0 to id
{
    gsl_matrix_view A=gsl_matrix_submatrix(matrix,k1,k2,n1,n2);
    Matrix* submat=new Matrix(n1,n2);
    gsl_matrix_memcpy(submat->matrix,&A.matrix);
    return submat;
}



double
Matrix::vectorMultiply(int i,double u1, int j,double u2)
{
    gsl_vector * x1 = gsl_vector_alloc(col);
    gsl_vector * x2 = gsl_vector_alloc(col);
    gsl_matrix_get_row(x1, matrix, i);
    gsl_matrix_get_row(x2, matrix, j);
    gsl_vector_add_constant (x1, -u1);
    //gsl_vector_view x2 = gsl_matrix_row(matrix, j);
    gsl_vector_add_constant (x2, -u2);
    double ssd1=0;
    gsl_blas_ddot(x1,x2,&ssd1);
    gsl_vector_free(x1);
    gsl_vector_free(x2);
    return ssd1;
}

// faster than gsl_stats_mean
double
Matrix::RowMean(int i)
{
    gsl_vector * x1 = gsl_vector_alloc(col);
    gsl_vector * x2 = gsl_vector_alloc(col);
    gsl_matrix_get_row(x1, matrix, i);
    gsl_vector_set_all(x2, 1.0);
    double ssd1=0;
    gsl_blas_ddot(x1,x2,&ssd1);
    gsl_vector_free(x1);
    gsl_vector_free(x2);
    double mean=ssd1/((double)col);
    return mean;
}

/*
double
Matrix::RowMean(int i)
{
    gsl_vector_view x1 = gsl_matrix_row(matrix, i);
    double mean=gsl_stats_mean(x1.vector.data,1,col);
    return mean;
}
 
//slow
double
Matrix::Covariance(int i,int j)
{
    gsl_vector_view x1 = gsl_matrix_row(matrix, i);
    gsl_vector_view x2 = gsl_matrix_row(matrix, j);
    double cov=gsl_stats_covariance(x1.vector.data, x1.vector.stride, x2.vector.data, x2.vector.stride, x1.vector.size);
    return cov;
}
 
//fast
double
Matrix::CovarianceGivMean(int i,int j,double imean,double jmean)
{
    gsl_vector_view x1 = gsl_matrix_row(matrix, i);
    gsl_vector_view x2 = gsl_matrix_row(matrix, j);
    double cov=gsl_stats_covariance_m(x1.vector.data, x1.vector.stride, x2.vector.data, x2.vector.stride, x1.vector.size, imean, jmean);
    return cov;
}*/




Matrix* 
Matrix::invMatrix() 		
{
	Matrix* minv=new Matrix(row,col);

	gsl_matrix* ludecomp=gsl_matrix_alloc(row,col);
	gsl_matrix_memcpy(ludecomp, matrix);
	//cout << "Old Value : " << gsl_matrix_get(matrix,0,0) << endl;
	//cout << "New Value : " << gsl_matrix_get(ludecomp,0,0) << endl;
	gsl_permutation* p=gsl_permutation_alloc(row);
	int signum=0;

	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	gsl_linalg_LU_invert(ludecomp, p, minv->matrix);
	gsl_matrix_free(ludecomp);
	gsl_permutation_free(p);
	return minv;
}

Matrix*
Matrix::invMatrix(gsl_matrix* ludecomp, gsl_permutation* p) 		
{
	//cout << "Old Value : " << gsl_matrix_get(matrix,0,0) << endl;
	//cout << "New Value : " << gsl_matrix_get(ludecomp,0,0) << endl;
	Matrix* minv=new Matrix(row,col);
	ludecomp->size1=row;
	ludecomp->size2=col;
	p->size=row;
	gsl_matrix_memcpy(ludecomp, matrix);
	int signum=0;
	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	gsl_linalg_LU_invert(ludecomp, p, minv->matrix);
	return minv;
}



Matrix*
Matrix::transMatrix()
{
	Matrix* transMatrix=new Matrix(col,row);
	gsl_matrix_transpose_memcpy(transMatrix->matrix,matrix);
	return transMatrix;
}

bool 
Matrix::dimequal(Matrix* aMatrix)
{
	if((aMatrix->getRowCnt()==row)&&
	   (aMatrix->getColCnt()==col))
	{
		return true;
	}	
	return false;
}

int 
Matrix::getRowCnt()
{
	return row;
}

int 
Matrix::getColCnt()
{
	return col;
}

double
Matrix::detMatrix()
{
	gsl_matrix* ludecomp=gsl_matrix_alloc(row,col);
	gsl_matrix_memcpy(ludecomp, matrix);
	gsl_permutation* p=gsl_permutation_alloc(row);
	int signum=0;

	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	double det=gsl_linalg_LU_det(ludecomp, signum);
	gsl_matrix_free(ludecomp);
	gsl_permutation_free(p);
	return det;
}

double
Matrix::detMatrix(gsl_matrix* ludecomp, gsl_permutation* p)
{
	ludecomp->size1=row;
	ludecomp->size2=col;
	p->size=row;
	gsl_matrix_memcpy(ludecomp, matrix);
	int signum=0;

	gsl_linalg_LU_decomp(ludecomp, p, &signum);
	double det=gsl_linalg_LU_det(ludecomp, signum);
	return det;
}



int 
Matrix::convertToFloat(gsl_matrix_float* dest,const gsl_matrix* source, int r, int c)
{
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			float val=gsl_matrix_get(source,i,j);	
			gsl_matrix_float_set(dest,i,j,val);
		}
	}
	return 0;
}

int 
Matrix::convertFromFloat(const gsl_matrix_float* source, gsl_matrix* dest,int r, int c)
{
	for(int i=0;i<r;i++)
	{
		for(int j=0;j<c;j++)
		{
			float val=gsl_matrix_float_get(source,i,j);
			gsl_matrix_set(dest,i,j,val);
		}
	}
	return 0;
}

int
Matrix::showMatrix(ostream& o)
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			o << getValue(i,j) << " ";
		}	
		o << endl;
	}
	return 0;
}

//Normalize column wise
int 
Matrix::normalize()
{
	double total=0;
	for(int i=0;i<col;i++)
	{
		for(int j=0;j<row;j++)
		{
			double val=getValue(j,i);
			total=total + val;
		}
	}	
	
	for(int i=0;i<col;i++)
	{
		for(int j=0;j<row;j++)
		{
			double val=getValue(j,i);
			val=val/total;	
			setValue(val,j,i);
		}	
	}
	return 0;
}

int
Matrix::normalizeVector()
{
	if(col>1)
	{
		return 0;
	}	
	double minValue=gsl_matrix_min(matrix);
	if(minValue<0)
	{
		minValue=-1*minValue;
		for(int i=0;i<row;i++)
		{
			double val=getValue(i,0);
			if(val<0)
			{
				val=val+minValue;
				setValue(val,i,0);	
			}
		}
	}

	double total=0;
	for(int i=0;i<row;i++)
	{
		total=total+getValue(i,0);
	}

	for(int i=0;i<row;i++)
	{
		double val=getValue(i,0)/total;
		setValue(val,i,0);
	}

	return 0;
}

double
Matrix::getMax()
{
	double val=gsl_matrix_max(matrix);
	return val;
}

int
Matrix::makeUncorrelated()
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(i!=j)
			{
				setValue(0,i,j);	
				//setValue(getValue(i,j),i,j);	
				//setValue(getValue(i,j),j,i);	
			}
			else 
			{
				if(getValue(i,j)==0)
				{
					cout << "What the f**k!! Zero variance u duffer!! at " << i << " : " << j <<endl; 
					setValue(0.001,i,j);
				}
			}
		}
	}
	return 0;
}

bool 
Matrix::rowZero()
{
	bool zeroCheck=false;
	for(int i=0;i<row;i++)
	{
		int zeroCol=0;
		for(int j=0;j<col;j++)
		{
			if(getValue(i,j)==0)	
			{
				zeroCol++;
			}
		}
		if(zeroCol==col)
		{
			cout << "Row " << i << " is zero" << endl;
			zeroCheck=true;
		}
	}
	return zeroCheck;
}

bool 
Matrix::colZero()
{
	bool zeroCheck=false;
	for(int i=0;i<col;i++)
	{
		int zeroRow=0;
		for(int j=0;j<row;j++)
		{
			if(getValue(j,i)==0)	
			{
				zeroRow++;
			}
		}
		if(zeroRow==row)
		{
			cout << "Col " << i << " is zero" <<  endl;
			zeroCheck=true;
		}
	}
	return zeroCheck;
}

//Make matrix non negative by making all elements less than zero
//as zero
int 
Matrix::makePositive()
{
	double minValue=gsl_matrix_min(matrix);
	if(minValue>=0)
	{
		//This means matrix is positive
		return 0;
	}

	minValue=minValue*-1;

	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			double value=getValue(i,j);
			value=value+minValue;
			setValue(value,i,j);
		}
	}
	return 0;	
}

Matrix*
Matrix::copyMe()
{
	Matrix* aMatrix=new Matrix(row,col);
	for(int i=0;i<row;i++)
	{	
		for(int j=0;j<col;j++)
		{
			double val=getValue(i,j);
			aMatrix->setValue(val,i,j);
		}
	}

	return aMatrix;
}

//Finds the closest vector to this vector
//which satisfies the constraint that it is a probability
//vector and also belongs to the same distribution

Matrix*
Matrix::findClosest()
{
	if(col!=1)
	{
		cout << "This is valid only for a vector " << endl;
		return NULL;
	}
	
	bool foundClosest=false;
	int iter=0;
	//Matrix* randVector=copyMe();
	Matrix* randVector=new Matrix(row,1);
	for(int i=0;i<row;i++)
	{
		double val=(double)rand()/RAND_MAX;
		randVector->setValue(val,i,0);
	}
	
	randVector->normalizeVector();
	double dist=0;

	while(iter < MAXITER)
	{
		//getDistance
		dist=getDistance(randVector);
		if(iter==0)
		{
			cout << "Initial Distance : " << dist << endl;
		}
		Matrix* nextClosest=getNextClosest(randVector,dist);
		if(nextClosest!=NULL)
		{
			delete randVector;
			randVector=nextClosest;
		}
		iter++;
	}	
	cout << "Final Distance : " << dist << endl;
	return randVector;
}

Matrix*
Matrix::getNextClosest(Matrix* current,double currDist)
{
	//Generate a vector for a list of possible closest
	//vectors. This is done by adding or subtracting a 
	//small fraction from the current 
	map<int,Matrix*> candidates;
	int id=0;
	double corr=(double)rand()/RAND_MAX;
	if(corr>0.1)
	{
		corr=corr*0.1;
	}
	
	for(int i=0;i<row;i++)
	{
		Matrix* cand1=current->copyMe();
		double val=cand1->getValue(i,0);
		double newVal=val+(val*corr);
		cand1->setValue(newVal,i,0);
		cand1->normalizeVector();
		Matrix* cand2=current->copyMe();
		newVal=val-(val*corr);
		cand2->setValue(newVal,i,0);
		cand2->normalizeVector();
		candidates[id]=cand1;
		id++;
		candidates[id]=cand2;
		id++;
	}

	double minDist=10000;
	int minVec=-1;
	for(int i=0;i<id;i++)
	{
		double dist=getDistance(candidates[i]);	
		if(dist<minDist)
		{
			minDist=dist;
			minVec=i;
		}
	}
	Matrix* bestFound=NULL;
	if(minDist<currDist)
	{
		bestFound=candidates[minVec];
	}
	
	map<int,Matrix*>::iterator anIter;
	for(anIter=candidates.begin();anIter!=candidates.end();anIter++)
	{
		Matrix* temp=anIter->second;
		candidates.erase(anIter);
		if(temp!=bestFound)
		{
			delete temp;
		}
	}
	
	return bestFound;
	
}

double
Matrix::getDistance(Matrix* a)
{
	double dist=0;
	for(int i=0;i<row;i++)
	{
		double val=a->getValue(i,0)-getValue(i,0);
		dist=dist + (val*val);
	}

	return sqrt(dist);
}
