LFLAG = -lgsl -lgslcblas 
SRC = Framework.C MetaLearner.C MetaMove.C SpeciesDataManager.C SpeciesDistance.C common/Error.C common/SlimFactor.C common/VariableManager.C common/Matrix.C common/PotentialManager.C common/Utils.C common/Variable.C

INCLPATH = common/

CC=g++
CFLAGS = -g -std=gnu++14 -gdwarf-2

condSpecLearner: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH) $(LFLAG) $(CFLAGS) -o scMTNI

