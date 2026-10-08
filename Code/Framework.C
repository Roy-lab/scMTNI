/*
 *  scMTNI: single-cell Multi-Task learning Network Inference
 *   Copyright 2022 Shilu Zhang (szhang256@wisc.edu) and  Sushmita Roy (sroy@biostat.wisc.edu)
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, 
 *   including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 *   subject to the following conditions:
 *
 *   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 *   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *   */
#include <fstream>
#include <iostream>
#include <cstring>

#include <unistd.h>

using namespace std;
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"
#include "MetaMove.H"
#include "SpeciesDistance.H"
#include "SpeciesDataManager.H"
#include "MetaLearner.H"
#include "Framework.H"

Framework::Framework()
{
    epsThreshold=-1;
}

Framework::~Framework()
{
}

Error::ErrorCode
Framework::init(int argc, char** argv)
{
    int optret='-';
    opterr=1;
    int oldoptind=optind;
    while(optret=getopt(argc,argv,"f:x:p:t:l:i:d:n:b:q:")!=-1)
    {
        if(optret=='?')
        {
            cout <<"Option error " << optopt << endl;
            return Error::UNKNOWN;
        }
        char c;
        char* my_optarg=NULL;
        c=*(argv[oldoptind]+1);
        if(optind-oldoptind ==2)
        {
            my_optarg=argv[oldoptind+1];
        }
        else
        {
            my_optarg=argv[oldoptind]+2;
        }
        switch(c)
        {
            case 'f':
            {
                metaLearner.setInputFName(my_optarg);
                break;
            }
            case 'x':
            {
                int aSize=atoi(my_optarg);
                metaLearner.setMaxFactorSize(aSize);
                break;
            }
            case 'p':
            {
                double proot=atof(my_optarg);
                spData.setProot(proot);
                break;
            }
            case 't':
            {
                double convThreshold=atof(my_optarg);
                metaLearner.setConvergenceThreshold(convThreshold);
                break;
            }
            case 'l':
            {
                metaLearner.setRegulatorList(my_optarg);
                break;
            }
            case 'i':
            {
                if(strcmp(my_optarg,"yes")==0)
                {
                    metaLearner.setINDEP();
                }
                break;
            }
            case 'd':
            {
                spData.readSpeciesTree(my_optarg);
                break;
            }
            case 'n':
            {
                metaLearner.setTargetList(my_optarg);
                break;
            }
            case 'b':
            {
                double beta1 = atof(my_optarg);
                metaLearner.setBeta1(beta1);
                break;
            }
            case 'q':
            {
              double beta2 = atof(my_optarg);
              metaLearner.setBeta2(beta2);
              break;
            }
            default:
            {
                cout <<"Unhandled option " << c  << endl;
                return Error::UNKNOWN;
            }
        }
        oldoptind=optind;
    }

    metaLearner.init();
    metaLearner.setSpeciesDistances(&spData);
    return Error::SUCCESS;
}

int 
Framework::start()
{
    metaLearner.start();
    return 0;
}

int
main(int argc, char* argv[])
{
    if(argc < 2) {
        cout <<"factorGraphInf " <<  endl
        << "-f modelname " << endl
        << "-x maxfactorsize" << endl
        << "-p proot" << endl
        << "-t convergence_threshold" << endl
        << "-l restrictedfname" << endl
        << "-i indep" << endl
        << "-d tree" << endl
        << "-n input_og_list" << endl
        << "-b beta1" << endl
        << "-q beta2" << endl;
        return 0;
    }

    Framework fw;
    if(fw.init(argc,argv)!=Error::SUCCESS) {
        return 0;
    }

    fw.start();
    return 0;
}

