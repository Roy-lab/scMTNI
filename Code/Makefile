LFLAG = -lgsl -lgslcblas 
SRC = Framework.C  MetaLearner.C  MetaMove.C	SpeciesDataManager.C	SpeciesDistance.C	 common/Error.C	   common/EvidenceManager.C  common/Potential.C	   common/SlimFactor.C	common/VariableManager.C common/Evidence.C  common/FactorGraph.C      common/Graph.C	     common/Matrix.C		common/PotentialManager.C	common/Utils.C  common/Variable.C	common/Vertex.C  orthomapsparser/MappedOrthogroupReader.C	orthomapsparser/MappedOrthogroup.C orthomapsparser/GeneMap.C
INCLPATH2 = common/
INCLPATH3 = orthomapsparser/

CC=g++
CFLAGS = -g -std=gnu++14 -gdwarf-2

condSpecLearner: $(SRC)
	$(CC) $(SRC) -I $(INCLPATH2) -I $(INCLPATH3)  $(LFLAG) $(CFLAGS) -o scMTNI



