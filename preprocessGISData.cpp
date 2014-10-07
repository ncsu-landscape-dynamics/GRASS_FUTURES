/****************************************************************************
 *
 * MODULE:       r.futures
 * AUTHOR(S):
 *
 * PURPOSE:
 *
 * COPYRIGHT:    (C) 2013-2014 by
 *
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <iostream>
#include <fstream>

extern "C" {
#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
}

#include "cfgutils.h"
//#include "distance.h"

#define	_CELL_OUT_OF_RANGE			-2			/* used to flag cells that are off the lattice */
#define	_CELL_OUT_OF_COUNTY			-1			/* used to flag cells that are not in this county */
#define	_CELL_VALID					1			/* used to flag cells that are valid */

#define _N_MAX_DYNAMIC_BUFF_LEN		(1024*1024)	/* need to dynamically allocate this as putting all on stack will crash most compilers */
#define _N_MAX_FILENAME_LEN			1024

#define	_GIS_HEADER_LENGTH			6			/* for working out where to start reading data */

#define	_GIS_NO_DATA_STRING			"-9999"		/* strictly-speaking, should probably parse these from GIS files */
#define	_GIS_NO_DATA_INT			-9999

#define	_N_NOT_YET_DEVELOPED		-1			/* use this for tDeveloped if cell still undeveloped */

/* algorithm to use */
#define	_N_ALGORITHM_DETERMINISTIC	1			/* deterministic model */
#define _N_ALGORITHM_STOCHASTIC_I	2			/* stochastic model with downweighting of logit probabilities and of devPressure over time */
#define	_N_ALGORITHM_STOCHASTIC_II	3

#define	_MAX_RAND_FIND_SEED_FACTOR	25
#define maxNumAddVariables          6
#define MAXNUM_COUNTY               50          /*maximal number of counties allowed*/
#define MAX_YEARS                   100
//#define MAX_UNDEV_SIZE              2000000     // maximum array size for undev cells (maximum: 1840198 for a county within 16 counties)
#define MAX_UNDEV_SIZE              2000000     // maximum array size for undev cells (maximum: 1840198 for a county within 16 counties)
using namespace std;


/* Wenwu Tang*/
char dirName[100]; // use absolute paths
char tempStr[100];// string for temporarily storage



typedef struct
{
	int			nCellType;						/* see #define's starting _CELL_ above */
	double		employAttraction;				/* attraction to employment base; static */
	double		interchangeDistance;			/* distance to interchange; static */
	double		roadDensity;					/* road density; static */
	int			thisX;							/* x position on the lattice; static */
	int			thisY;							/* y position on the lattice; static */
	int			bUndeveloped;					/* whether this site is still undeveloped; varies.  Note if bUndeveloped = 0 the values of employAttraction etc. are not to be trusted */
	double		devPressure;					/* development pressure; varies */
	int			bUntouched;						/* stores whether this site has been considered for development, to handle "protection" */
	int			tDeveloped;						/* timestep on which developed (0 for developed at start, _N_NOT_YET_DEVELOPED for not developed yet ) */
	double		consWeight;						/* multiplicative factor on the probabilities */
    double      additionVariable[maxNumAddVariables];               /* additional variables*/
    int         index_region;
    float       devProba;
} t_Cell;

typedef struct
{
	int			cellID;							/* id of this cell */
	double		logitVal;						/* value of the logit */
	int			bUntouched;						/* whether or not cell previously considered...need to use this in the sort */
	double      cumulProb;                      /* to support random pick based on their logitVal*/
} t_Undev;

typedef struct
{
	t_Cell		*asCells;						/* array of cells (see posFromXY for how they are indexed) */
	int			maxX;							/* number of columns in the grid */
	int			maxY;							/* number of rows in the grid */
	int			totalCells;						/* total number of cells */
	t_Undev		*asUndev;						/* array of information on undeveloped cells */
	int			undevSites;						/* number of cells which have not yet converted */
	int			*aParcelSizes;					/* posterior sample from parcel size distribution */
	int			parcelSizes;					/* number in that sample */
    t_Undev     **asUndevs;                     //WT
    int         num_undevSites[MAXNUM_COUNTY];  //WT
} t_Landscape;

typedef struct
{
	/* size of the grids (could be worked out, but easier this way) */
	int			xSize;											/* 7312 */
	int			ySize;											/* 5571 */
	/* file containing information on how many cells to transition and when */
    char*		controlFile;
	/* files containing the information to read in */
    char*		employAttractionFile;		/* atr1500clip */
    char*		interchangeDistanceFile;	/* d2interclip */
    char*		roadDensityFile;			/* psden1000clip */
    char*		undevelopedFile;			/* cab_undv_msk */
    char*		devPressureFile;			/* devp500clip */
    char*		consWeightFile;			/* consweight */
    char*		probLookupFile;
	int			nProbLookup;
	double		*adProbLookup;
	/* parameters that go into regression formula */
	double		dIntercept;										/* 0.038884 */
	double		dEmployAttraction;								/* -0.0000091946 */
	double		dInterchangeDistance;							/* 0.000042021 */
	double		dRoadDensity;									/* 0.00065813 */
	double		dDevPressure;									/* -0.026190 */
	/* size of square used to recalculate development pressure */
	int			nDevNeighbourhood;								/* 8 (i.e. 8 in each direction, leading to 17 x 17 = 289 in total) */
	/* used in writing rasters */
    char*		dumpFile;
	/* 1 deterministic, 2 old stochastic, 3 new stochastic */
	int			nAlgorithm;
	/* use this to downweight probabilities */
	double		dProbWeight;
	/* and this to downweight old dev pressure */
	double		dDevPersistence;
	/* file containing parcel size information */
    char*		parcelSizeFile;
       double        discountFactor;//for calibrate patch size
	/* give up spiralling around when examined this number of times too many cells */
	double		giveUpRatio;
	/* these parameters only relevant for new stochastic algorithm */
	int			sortProbs;
	double		patchFactor;
	double      patchMean;
	double      patchRange;
	int         numNeighbors; // 4 or 8 neighbors only
	int         seedSearch; //1: uniform distribution 2: based on dev. proba.
	int         numAddVariables;
	double      addParameters[maxNumAddVariables][MAXNUM_COUNTY]; /* parameters for additional variables*/
    char *       addVariableFile[maxNumAddVariables];
    int         devPressureApproach; //1: #occurrence; 2: gravity (need to define alpha and scaling factor); 3: kernel, need to define alpha and scaling factor
                                     // for 2: formula-> scalingFactor/power(dist,alpha)
                                     // for 3: formula-> scalingFactor*exp(-2*dist/alpha)
    double      scalingFactor; // scale distance-based force
    double      alpha; // constraint on distance

    /* parameters that go into regression formula */
	double		dIntercepts[MAXNUM_COUNTY];										/* 0.038884 */
	double		dV1[MAXNUM_COUNTY];								/* -0.0000091946 */
	double		dV2[MAXNUM_COUNTY];							/* 0.000042021 */
	double		dV3[MAXNUM_COUNTY];									/* 0.00065813 */
	double		dV4[MAXNUM_COUNTY];
    int         num_Regions;
    char		indexFile[_N_MAX_FILENAME_LEN];	/* index file to run multiple regions */
    /*control files for all regions*/
    char		controlFileAll[_N_MAX_FILENAME_LEN]; //development demand for multiple regions
    int         devDemand[MAX_YEARS];
    int         devDemands[MAXNUM_COUNTY][MAX_YEARS];
    int         nSteps; //#simulation steps
} t_Params;


//From WT: move here to be global variables
t_Params	sParams;
t_Landscape	sLandscape;

int getUnDevIndex(t_Landscape *pLandscape);
void Shuffle4(int*mat);
void print2ASC(t_Landscape*pLandscape,char* fn);
double getDistance1(double x1,double y1,double x2,double y2);
void readDevPotParams(t_Params *pParams,char*fn);
void readIndexData(t_Landscape* pLandscape, t_Params *pParams);
void findAndSortProbsAll(t_Landscape *pLandscape, t_Params *pParams,int step);
void updateMap1(t_Landscape *pLandscape, t_Params *pParams, int step, int regionID);
int getUnDevIndex1(t_Landscape *pLandscape,int regionID);
void export2ASC1(t_Landscape *pLandscape, t_Params *pParams, int regionID,char* fn);

char* addPath(char*str1,char* pathStr){
	//add directory name directly here to use absolute path // by Wenwu Tang
    char tempS[100];
    strcpy(tempS,str1);
    strcpy(str1,pathStr);
    strcat(str1,tempS);
    return str1;
}

/*
	read in parameters from specified config file
*/
int setParams(t_Params *pParams, char *szCfgFile)
{
	fprintf(stdout, "entered readParams()\n");


	/* size of the grids (could be worked out, but easier this way) */
	if(!readIntFromCfg(szCfgFile, "xSize", &pParams->xSize))
	{
		fprintf(stderr, "error reading xSize...exiting\n");
		return 0;
	}
	if(!readIntFromCfg(szCfgFile, "ySize", &pParams->ySize))
	{
		fprintf(stderr, "error reading ySize...exiting\n");
		return 0;
	}

	/* file containing information on how many cells to transition and when */
	if(!readStringFromCfg(szCfgFile, "controlFile", pParams->controlFile))
	{
		fprintf(stderr, "error reading controlFile...exiting\n");
		return 0;
	}
    addPath(pParams->controlFile,dirName);
	/* files containing the information to read in */
	if(!readStringFromCfg(szCfgFile, "employAttractionFile", pParams->employAttractionFile))
	{
		fprintf(stderr, "error reading employAttractionFile...exiting\n");
		return 0;
	}
	addPath(pParams->employAttractionFile,dirName);
	if(!readStringFromCfg(szCfgFile, "interchangeDistanceFile", pParams->interchangeDistanceFile))
	{
		fprintf(stderr, "error reading interchangeDistanceFile...exiting\n");
		return 0;
	}
	addPath(pParams->interchangeDistanceFile,dirName);
	if(!readStringFromCfg(szCfgFile, "roadDensityFile", pParams->roadDensityFile))
	{
		fprintf(stderr, "error reading roadDensityFile...exiting\n");
		return 0;
	}
	addPath(pParams->roadDensityFile,dirName);
	if(!readStringFromCfg(szCfgFile, "undevelopedFile", pParams->undevelopedFile))
	{
		fprintf(stderr, "error reading undevelopedFile...exiting\n");
		return 0;
	}
	addPath(pParams->undevelopedFile,dirName);
	if(!readStringFromCfg(szCfgFile, "devPressureFile", pParams->devPressureFile))
	{
		fprintf(stderr, "error reading devPressureFile...exiting\n");
		return 0;
	}
	addPath(pParams->devPressureFile,dirName);
	if(!readStringFromCfg(szCfgFile, "consWeightFile", pParams->consWeightFile))
	{
		fprintf(stderr, "error reading consWeightFile...exiting\n");
		return 0;
	}
	addPath(pParams->consWeightFile,dirName);

	/* parameters that go into regression formula */
    /* not used
	if(!readDoubleFromCfg(szCfgFile, "dIntercept", &pParams->dIntercept))
	if(!readDoubleFromCfg(szCfgFile, "dEmployAttraction", &pParams->dEmployAttraction))
	if(!readDoubleFromCfg(szCfgFile, "dInterchangeDistance", &pParams->dInterchangeDistance))
	if(!readDoubleFromCfg(szCfgFile, "dRoadDensity", &pParams->dRoadDensity))
	if(!readDoubleFromCfg(szCfgFile, "dDevPressure", &pParams->dDevPressure))
    */
    /*Wenwu add additional variables*/
    if(!readIntFromCfg(szCfgFile, "numAddVariables", &pParams->numAddVariables))
	{
		fprintf(stderr, "error reading numAdditionalVariables...exiting\n");
		return 0;
	}
	char str[50],tempStr[20]; int i;
    for(i=0;i<pParams->numAddVariables;i++){
        //read parameters
        strcpy(str,"parameterLogistic");
        sprintf(tempStr,"%d",i+1);
        strcat(str,tempStr);
        /*
        if(!readDoubleFromCfg(szCfgFile, str, &pParams->addParameters[i]))
        {
            fprintf(stderr, "error reading additionalParameters...exiting\n");
            return 0;
        }
        */
        //read file names
        strcpy(str,"variableFile");
        strcat(str,tempStr);
        if(!readStringFromCfg(szCfgFile, str, pParams->addVariableFile[i]))
        {
            fprintf(stderr, "error reading additionalVariableFile...exiting\n");
            return 0;
        }
        addPath(pParams->addVariableFile[i],dirName);
    }

	/* size of square used to recalculate development pressure */
	if(!readIntFromCfg(szCfgFile, "nDevNeighbourhood", &pParams->nDevNeighbourhood))
	{
		fprintf(stderr, "error reading nDevNeighbourhood...exiting\n");
		return 0;
	}

	/* used in writing rasters */
	if(!readStringFromCfg(szCfgFile, "dumpFile", pParams->dumpFile))
	{
		fprintf(stderr, "error reading dumpFile...exiting\n");
		return 0;
	}
    addPath(pParams->dumpFile,dirName);
	/* parameters controlling the algorithm to use */
    if(!readIntFromCfg(szCfgFile, "nAlgorithm", &pParams->nAlgorithm))
	{
		fprintf(stderr, "error reading nAlgorithm...exiting\n");
		return 0;
	}
    if(!readDoubleFromCfg(szCfgFile, "dProbWeight", &pParams->dProbWeight))
	{
		fprintf(stderr, "error reading dProbWeight...exiting\n");
		return 0;
	}
	if(!readDoubleFromCfg(szCfgFile, "dDevPersistence", &pParams->dDevPersistence))
	{
		fprintf(stderr, "error reading dDevPersistence...exiting\n");
		return 0;
	}

	/* file containing information on the parcel size to use */
	if(!readStringFromCfg(szCfgFile, "parcelSizeFile", pParams->parcelSizeFile))
	{
		fprintf(stderr, "error reading parcelSizeFile...exiting\n");
		return 0;
	}
	addPath(pParams->parcelSizeFile,dirName);
       if(!readDoubleFromCfg(szCfgFile, "discountFactor", &pParams->discountFactor))
	{
		fprintf(stderr, "error reading discountFactor of patch size...exiting\n");
		return 0;
	}    
	if(!readDoubleFromCfg(szCfgFile, "giveUpRatio", &pParams->giveUpRatio))
	{
		fprintf(stderr, "error reading giveUpRatio...exiting\n");
		return 0;
	}
	pParams->sortProbs = 1;
	if(pParams->nAlgorithm == _N_ALGORITHM_STOCHASTIC_II)
	{
		int	 parsedOK,i;
		FILE *fp;
		char inBuff[N_MAXREADINLEN];
		char *pPtr;

		fprintf(stdout, "reading probability lookup\n");
		if(!readStringFromCfg(szCfgFile, "probLookup", pParams->probLookupFile))
		{
			fprintf(stderr, "error reading probLookup...exiting\n");
			return 0;
		}
        addPath(pParams->probLookupFile,dirName);

		fp = fopen(pParams->probLookupFile,"rb");
		if(fp)
		{
			parsedOK = 0;
			if(fgets(inBuff,N_MAXREADINLEN,fp))
			{
				if(inBuff[0] == ',')
				{
					pParams->nProbLookup = atoi(inBuff+1);
					if(pParams->nProbLookup > 0)
					{
						pParams->adProbLookup = (double*)malloc(sizeof(double)*pParams->nProbLookup);
						if(pParams->adProbLookup)
						{
							parsedOK = 1;
							i=0;
							while(parsedOK && i < pParams->nProbLookup)
							{
								parsedOK = 0;
								if(fgets(inBuff,N_MAXREADINLEN,fp))
								{
									if(pPtr = strchr(inBuff, ','))
									{
										parsedOK=1;
										pParams->adProbLookup[i] = atof(pPtr + 1);
/*										fprintf(stdout, "\t%d %f->%f\n", i, (double)i*1.0/(pParams->nProbLookup-1), pParams->adProbLookup[i]); */
									}
								}
								i++;
							}
						}
					}
				}
			}
			if(!parsedOK)
			{
				fprintf(stderr, "error parsing probLookup file '%s'...exiting\n", pParams->probLookupFile);
				return 0;
			}
			fclose(fp);
		}
		else
		{
			fprintf(stderr, "error opening probLookup file '%s'...exiting\n", pParams->probLookupFile);
			return 0;
		}
		/* parameters controlling the newer algorithm */
		if(!readIntFromCfg(szCfgFile, "sortProbs", &pParams->sortProbs))
		{
			fprintf(stderr, "error reading sortProbs...exiting\n");
			return 0;
		}
		if(!readDoubleFromCfg(szCfgFile, "patchFactor", &pParams->patchFactor))
		{
			fprintf(stderr, "error reading patchFactor...exiting\n");
			return 0;
		}
		if(!readDoubleFromCfg(szCfgFile, "patchMean", &pParams->patchMean))
		{
			fprintf(stderr, "error reading patchMean...exiting\n");
			return 0;
		}
		if(!readDoubleFromCfg(szCfgFile, "patchRange", &pParams->patchRange))
		{
			fprintf(stderr, "error reading patchRange...exiting\n");
			return 0;
		}
        if(!readIntFromCfg(szCfgFile, "numNeighbors", &pParams->numNeighbors))
		{
			fprintf(stderr, "error reading numNeighbors...exiting\n");
			return 0;
		}
        if(!readIntFromCfg(szCfgFile, "seedSearch", &pParams->seedSearch))
		{
			fprintf(stderr, "error reading seedSearch...exiting\n");
			return 0;
		}
        if(!readIntFromCfg(szCfgFile, "devPressureApproach", &pParams->devPressureApproach))
        {
            fprintf(stderr, "error reading devPressureApproach...exiting\n");
            return 0;
        }
        if(pParams->devPressureApproach!=1){
        if(!readDoubleFromCfg(szCfgFile, "alpha", &pParams->alpha))
        {
            fprintf(stderr, "error reading alpha...exiting\n");
            return 0;
        }
        if(!readDoubleFromCfg(szCfgFile, "scalingFactor", &pParams->scalingFactor))
        {
            fprintf(stderr, "error reading scalingFactor...exiting\n");
            return 0;
        }
        }
        if(!readIntFromCfg(szCfgFile, "num_regions", &pParams->num_Regions))
        {
            fprintf(stderr, "error reading devPressureApproach...exiting\n");
            return 0;
        }
        if(pParams->num_Regions>1)
            readDevPotParams(pParams,"./devpotParams.cfg");
        if(!readStringFromCfg(szCfgFile, "indexFile", pParams->indexFile))
        {
            fprintf(stderr, "error reading indexFile...exiting\n");
            return 0;
        }
        addPath(pParams->indexFile,dirName);
        if(!readStringFromCfg(szCfgFile, "controlFileAll", pParams->controlFileAll))
        {
            fprintf(stderr, "error reading control File All...exiting\n");
            return 0;
        }
        addPath(pParams->controlFileAll,dirName);
	}
	return 1;
}
void readDevPotParams(t_Params *pParams,char*fn){
    char str[200]; int id; double di,d1,d2,d3,d4,d5,d6,val;
    int i,j;
    ifstream f;
    f.open(fn);
    f.getline(str,200);
    for(i=0;i<pParams->num_Regions;i++){
        f>>id>>di>>d1>>d2>>d3>>d4;
cout<<id<<"\t"<<di<<"\t"<<d1<<"\t"<<d2<<"\t"<<d3<<"\t"<<d4<<endl;

        pParams->dIntercepts[i]=di;
        pParams->dV1[i]=d1;
        pParams->dV2[i]=d2;
        pParams->dV3[i]=d3;
        pParams->dV4[i]=d4;
        for(j=0;j<pParams->numAddVariables;j++){
            f>>val;
            pParams->addParameters[j][i]=val;
        }
    }
    f.close();
}
/*
	return uniform number on [0,1)

	encapsulated to allow easy replacement if necessary
*/
double	uniformRandom()
{
	int nRet;

	nRet = RAND_MAX;
	while(nRet == RAND_MAX)			/* make sure never get 1.0 */
	{
		nRet = rand();
	}
	return ((double)nRet/(double)(RAND_MAX));
}

/*
	seed random number generator

	encapsulated to allow easy replacement if necessary
*/
void	seedRandom(struct timeval ttime)
{
	srand((ttime.tv_sec * 100) + (ttime.tv_usec / 100));
	//srand((unsigned int)time(NULL));
}


/*
	work out x and y position on the grid from index into siteMapping array
*/
void xyFromPos(int nPos, int *pnX, int *pnY, t_Landscape *pLandscape)
{
	*pnX = nPos % pLandscape->maxX;
	*pnY = nPos / pLandscape->maxX;	/* integer division just gives whole number part */
}

/*
	work out index into asCells array from location
*/
int	posFromXY(int nX, int nY, t_Landscape *pLandscape)
{
	int nPos;

	/*
		0   1	 2 		.... nX-1
		nX  nX+1 nX+2   .... 2*nX-1
		.						.
		.						.
		.						.
		(nY-1)*nX		...  nX*nY-1
	*/
	nPos = _CELL_OUT_OF_RANGE;
	if(nX >= 0 && nX < pLandscape->maxX && nY >= 0 && nY < pLandscape->maxY)
	{
		nPos = nX + nY * pLandscape->maxX;
	}
	return nPos;
}

/*
	allocate memory to store information on cells
*/
int buildLandscape(t_Landscape *pLandscape, t_Params *pParams)
{
	int bRet;

	/* blank everything out */
	fprintf(stdout, "entered buildLandscape()\n");
	bRet = 0;
	memset(pLandscape,0,sizeof(t_Landscape));

	/* set initial values across the landscape */
	pLandscape->undevSites = 0;
	pLandscape->maxX = pParams->xSize;
	pLandscape->maxY = pParams->ySize;
	pLandscape->totalCells = pLandscape->maxX * pLandscape->maxY;
	pLandscape->asCells = (t_Cell*)malloc(sizeof(t_Cell) * pLandscape->totalCells);
	if(pLandscape->asCells)
	{
	    if(pParams->num_Regions==1){
            pLandscape->asUndev = (t_Undev*)malloc(sizeof(t_Undev) * pLandscape->totalCells);	/* just allocate enough space for every cell to be undev */
            if(pLandscape->asUndev)
            {
                bRet = 1;
            }
        }
        else
            bRet=1;
	}
	return bRet;
}

void readData4AdditionalVariables(t_Landscape* pLandscape, t_Params *pParams){
    int i,j; char str[50];
    ifstream f;
    int ii,jj; double val;
    for(i=0;i<pParams->numAddVariables;i++){
        f.open(pParams->addVariableFile[i]);
        cout<<"reading additional variables File: "<<pParams->addVariableFile[i]<<"...";
        for(j=0;j<6;j++) f.getline(str,50);
        for(ii=0;ii<pParams->xSize*pParams->ySize;ii++){
            f>>val;
            if(pLandscape->asCells[ii].nCellType == _CELL_VALID){
                pLandscape->asCells[ii].additionVariable[i]=val;
            }
            else pLandscape->asCells[ii].additionVariable[i]=0.0;
        }
        f.close();
        cout<<"done"<<endl;
    }
}
void readIndexData(t_Landscape* pLandscape, t_Params *pParams){
    int i,j; char str[50];
    ifstream f;
    int ii,jj; int val;
    f.open(pParams->indexFile);
    cout<<"reading index File: "<<pParams->indexFile<<"...";
    for(j=0;j<6;j++) f.getline(str,50);
    for(ii=0;ii<pParams->xSize*pParams->ySize;ii++){
        f>>val;
        pLandscape->asCells[ii].index_region=val;
        if(val==3)
            int stop=1;

    }
    f.close();
    cout<<"done"<<endl;

}
/*
	unavoidably unpleasant routine to read data from GIS rasters and put in correct places in memory
*/
int	readData(t_Landscape *pLandscape, t_Params *pParams)
{
	FILE 	*fIn;
	char	*szBuff,*pPtr;
	char	szFName[_N_MAX_FILENAME_LEN];
	int		bRet,x,y,i,j,headerCount;
	double	dVal;

	fprintf(stdout, "entered readData()\n");
	bRet = 0;
	szBuff = (char*)malloc(_N_MAX_DYNAMIC_BUFF_LEN*sizeof(char));
	if(szBuff)
	{
		for(j=0;j<6;j++)
		{
			switch(j)	/* get correct filename */
			{
			case 0:
				strcpy(szFName, pParams->undevelopedFile);
				break;
			case 1:
				strcpy(szFName, pParams->employAttractionFile);
				break;
			case 2:
				strcpy(szFName, pParams->interchangeDistanceFile);
				break;
			case 3:
				strcpy(szFName, pParams->roadDensityFile);
				break;
			case 4:
				strcpy(szFName, pParams->devPressureFile);
				break;
			case 5:
				strcpy(szFName, pParams->consWeightFile);
				break;
			default:
				fprintf(stderr, "readData(): shouldn't get here...\n");
				break;
			}
			fprintf(stdout, "\t%s...", szFName);
			fIn = fopen(szFName,"rb");
			if(fIn)
			{
				/* strip GIS header */
				headerCount = 0;
				while(headerCount < _GIS_HEADER_LENGTH && fgets(szBuff, _N_MAX_DYNAMIC_BUFF_LEN, fIn))
				{
					headerCount++;
				}
				i = 0;
				y = 0;
				bRet = 1;	/* can only get worse from here on in */
				while(bRet && fgets(szBuff, _N_MAX_DYNAMIC_BUFF_LEN, fIn)) 	/* will fail on long lines */
				{
					if(y >=	pLandscape->maxY)
					{
						if(bRet)
						{
							fprintf(stderr, "readData(): y too large\n");
							bRet = 0;
						}
					}
					pPtr = strpbrk(szBuff, "\r\n");
					if(pPtr)
					{
						pPtr[0] = '\0';
					}
					else
					{
						fprintf(stderr, "readData(): line too long\n");
						bRet=0;				/* no CR or NL means line not read in fully */
					}
					pPtr = strtok(szBuff, " \t");
					x = 0;
					while(bRet && pPtr)
					{
						dVal = atof(pPtr);
						if(x < pLandscape->maxX)
						{
							if(j == 0)	/* check for NO_DATA */
							{
								if(strcmp(pPtr, _GIS_NO_DATA_STRING)==0)
								{
									pLandscape->asCells[i].nCellType = _CELL_OUT_OF_COUNTY;
								}
								else
								{
									pLandscape->asCells[i].nCellType = _CELL_VALID;
								}
								/* and put in the correct place */
								pLandscape->asCells[i].thisX = x;
								pLandscape->asCells[i].thisY = y;
								/* and set the time of development */
								pLandscape->asCells[i].tDeveloped = _GIS_NO_DATA_INT;
							}
							else
							{
								if(strcmp(pPtr, _GIS_NO_DATA_STRING)==0)	/* clean up missing data */
								{
									dVal = 0.0;
								}
							}
							if(pLandscape->asCells[i].nCellType == _CELL_VALID)
							{
								switch(j)	/* put in correct place */
								{
								case 0:
									pLandscape->asCells[i].bUndeveloped = atoi(pPtr);
									pLandscape->asCells[i].bUntouched = 1;
									if(pLandscape->asCells[i].bUndeveloped == 1)
									{
										pLandscape->asCells[i].tDeveloped = _N_NOT_YET_DEVELOPED;
									}
									else
									{
										pLandscape->asCells[i].tDeveloped = 0;	/* already developed at the start */
									}
									break;
								case 1:
									pLandscape->asCells[i].employAttraction = dVal;
									break;
								case 2:
									pLandscape->asCells[i].interchangeDistance = dVal;
									break;
								case 3:
									pLandscape->asCells[i].roadDensity = dVal;
									break;
								case 4:
									pLandscape->asCells[i].devPressure = (int) dVal;
									break;
								case 5:
									pLandscape->asCells[i].consWeight = dVal;
									break;
								default:
									fprintf(stderr, "readData(): shouldn't get here...\n");
									break;
								}
							}
						}
						else
						{
							bRet = 0;
							fprintf(stderr, "readData(): x too large\n"); /* too large in x direction */
						}
						i++;
						x++;
						pPtr = strtok(NULL, " \t");
					}
					if(x==pLandscape->maxX)
					{
						y++;
					}
					else
					{
						if(bRet)
						{
							fprintf(stderr, "readData(): x too small\n");
							bRet = 0;
						}
					}
				}
				if(y!=pLandscape->maxY)
				{
					if(bRet)
					{
						fprintf(stderr, "readData(): y too small\n");
						bRet = 0;
					}
				}
				else
				{
					if(bRet && i == pLandscape->totalCells)
					{
						fprintf(stdout, "done\n");;
					}
					else
					{
						if(bRet)
						{
							fprintf(stderr, "readData(): not read in enough points\n");
							bRet = 0;
						}
					}
				}
				fclose(fIn);
			}
			else
			{
				fprintf(stderr, "readData(): couldn't open %s\n", szFName);
			}
		}
		free(szBuff);
	}
	return bRet;
}


/*
	used in qsort...note that this sort function is in ascending order (i.e. opposite to what we normally want)
*/
static int undevCmp(const void *a, const void *b)
{
	double da = ((t_Undev *)a)->logitVal;
	double db = ((t_Undev *)b)->logitVal;

	/* make sure untouched cells preferentially go to top of sort when there is asymmetry between a and b */
	if(((t_Undev *)a)->bUntouched && !((t_Undev *)b)->bUntouched)
		return -1;
	if(((t_Undev *)b)->bUntouched && !((t_Undev *)a)->bUntouched)
		return 1;
	/* then sort by value */
	if(da > db)
		return 1;
	if(da < db)
		return -1;
	return 0;
}

/*
	note this is the correct way around...
*/
static int undevCmpReverse(const void *a, const void *b)
{
	return -1 * undevCmp(a,b);
}


/*
	called at end and dumps tDeveloped for all valid cells (-9999 for all others)
*/
void dumpDevRaster(t_Landscape *pLandscape, t_Params *pParams, char *szFile)
{
	FILE	*fIn,*fOut;
	char	*szBuff,*pPtr;
	int		x,y,cellID,toPrint,bVal;

	fIn = fopen(pParams->employAttractionFile, "rb");
	if(fIn)
	{
		fOut = fopen(szFile, "wb");
		if(fOut)
		{
			szBuff = (char*) malloc(_N_MAX_DYNAMIC_BUFF_LEN*sizeof(char));
			if(szBuff)
			{
				/* copy header across from one of the inputs */
				y=0;
				while(y<_GIS_HEADER_LENGTH && fgets(szBuff, _N_MAX_DYNAMIC_BUFF_LEN, fIn))
				{
					pPtr = strpbrk(szBuff, "\r\n");
					if(pPtr)
					{
						pPtr[0] = '\0';
					}
					fprintf(fOut, "%s\n", szBuff);
					y++;
				}
				/* now dump actual information */
				for(y=0;y<pLandscape->maxY;y++)
				{
					for(x=0;x<pLandscape->maxX;x++)
					{
						toPrint = _GIS_NO_DATA_INT;
						bVal = 0;
						cellID = posFromXY(x, y, pLandscape);
						if(cellID != _CELL_OUT_OF_RANGE)	/* should never happen */
						{
							if(pLandscape->asCells[cellID].nCellType == _CELL_VALID)
							{
								toPrint = pLandscape->asCells[cellID].tDeveloped;
								bVal = 1;
							}
						}
						if(x)
						{
							fprintf(fOut, " ");
						}
						if(bVal)
						{
							fprintf(fOut, "%d", toPrint);
						}
						else
						{
							fprintf(fOut, "%s", _GIS_NO_DATA_STRING);
						}
					}
					fprintf(fOut, "\n");
				}
				free(szBuff);
			}
			fclose(fOut);
		}
		fclose(fIn);
	}
}

/*
	work out how many cells to transition based on the control file
*/
int parseControlLine(char *szBuff,char *szStepLabel,int *pnToConvert)
{
	char *pPtr;
	int	 i;
	int bRet;

	bRet = 0;
	/* strip newline */
	pPtr = strpbrk(szBuff, "\r\n");
	if(pPtr)
	{
		pPtr[0] = '\0';
	}
	*pnToConvert = 0;
	i = 0;
	pPtr = strtok(szBuff, " \t");
	while(pPtr && i < 2)
	{
		switch(i)
		{
		case 0:
			strcpy(szStepLabel, pPtr);
			break;
		case 1:
			*pnToConvert = atoi(pPtr);
			break;
		}
		pPtr = strtok(NULL, " \t");
		i++;
	}
	if(i == 2)	/* only process this line of the control file if it is well formed */
	{
		bRet = 1;
	}
	return bRet;
}
/*
calculate probability
*/

double getDevProbability(t_Cell	*pThis, t_Params *pParams){
    float probAdd; int i;
    int id=pThis->index_region;
    if(id==-9999)
        return 0;
    id=id-1;
    probAdd = pParams->dIntercepts[id];
//cout<<"intercept\t"<<probAdd<<endl;
    probAdd += pParams->dV1[id] * pThis->employAttraction;
//cout<<"employAttraction: "<<pParams->dV1[id]<<endl;
//cout<<"employAttraction: "<<pThis->employAttraction<<endl;
//cout<<"employAttraction: "<<pParams->dV1[id] * pThis->employAttraction<<endl;
    probAdd += pParams->dV2[id] * pThis->interchangeDistance;
//cout<<"interchangeDistance: "<<pParams->dV2[id]<<endl;
//cout<<"interchangeDistance: "<<pThis->interchangeDistance<<endl;
//cout<<"interchangeDistance: "<<pParams->dV2[id] * pThis->interchangeDistance<<endl;
    probAdd += pParams->dV3[id] * pThis->roadDensity;
//cout<<"roadDensity: "<<pParams->dV3[id]<<endl;
//cout<<"roadDensity: "<<pThis->roadDensity<<endl;
//cout<<"roadDensity: "<<pParams->dV3[id] * pThis->roadDensity<<endl;
    probAdd += pParams->dV4[id] * pThis->devPressure;
//cout<<"devPressure: "<<pParams->dV4[id]<<endl;
//cout<<"devPressure: "<<pThis->devPressure<<endl;
//cout<<"devPressure: "<<pParams->dV4[id] * pThis->devPressure<<endl;
    for(i=0;i<pParams->numAddVariables;i++){
        probAdd += pParams->addParameters[i][id] * pThis->additionVariable[i];
//cout<<"additionVariable: "<<i<<"\t"<<pParams->addParameters[i][id]<<endl;
//cout<<"additionVariable: "<<i<<"\t"<<pThis->additionVariable[i]<<endl;
//cout<<"additionVariable: "<<i<<"\t"<<pParams->addParameters[i][id] * pThis->additionVariable[i]<<endl;
    }
	
    probAdd = 1.0 / (1.0 + exp(-probAdd));
//cout<<"The probability:";
//cout<<probAdd<<endl;
    return probAdd;
}

/*
double getDevProbability(t_Cell	*pThis, t_Params *pParams){
    float probAdd; int i;
    int id=pThis->index_region;
    if(id==-9999)
        return 0;
    id=id-1;
    probAdd = pParams->dIntercepts[id];
    probAdd += pParams->dV1[id] * pThis->employAttraction;
    probAdd += pParams->dV2[id] * pThis->interchangeDistance;
    probAdd += pParams->dV3[id] * pThis->roadDensity;
    probAdd += pParams->dV4[id] * pThis->devPressure;
    for(i=0;i<2;i++){
        probAdd+=pParams->addParameters[i][id]*pThis->additionVariable[i];
    }
	if (0<pThis->additionVariable[2]<=400)
	{
	probAdd = probAdd * 0.0005 ;
	}
	  else { 
	     if(400<pThis->additionVariable[2]<=600)
		 {probAdd = probAdd * 0.024;}
		     else {
			   if(600<pThis->additionVariable[2]<=800)
			        {probAdd = probAdd * 0.912 ;}
			       else{
				       if(800<pThis->additionVariable[2]<=1000)
					      {probAdd = probAdd * 0.055 ;}
						   else{
							if(1000<pThis->additionVariable[2]<=1200)
						    {probAdd = probAdd * 0.006 ;}
							    else{
								if(1200<pThis->additionVariable[2]<=1400)
								 {probAdd = probAdd * 0.002 ;}
								  else{
								   
									 {probAdd = probAdd * 0.0005;}

									 }
								  }
								}
						   
						   }
						   
				   }
			    
			 
	}
	
	
    probAdd = 1.0/(1.0 + exp(-probAdd));
    return probAdd;
}

*/


/*
	called each timestep...recalculate probabilities for each unconverted cell
*/
void findAndSortProbs(t_Landscape *pLandscape, t_Params *pParams, int nToConvert)
{
	int			i,lookupPos;
	t_Cell		*pThis;

	/* update calcs */
	fprintf(stdout, "\t\trecalculating probabilities\n");
	pLandscape->undevSites = 0;
	for(i=0;i<pLandscape->totalCells;i++)
	{
		pThis = &(pLandscape->asCells[i]);
		if(pThis->nCellType ==_CELL_VALID)
		{
			if(pThis->bUndeveloped)
			{
				if(pThis->consWeight > 0.0)
				{
					/* note that are no longer just storing the logit value, but instead the probability (allows consWeight to affect sort order) */
					pLandscape->asUndev[pLandscape->undevSites].cellID = i;
					pLandscape->asUndev[pLandscape->undevSites].logitVal=getDevProbability(pThis,pParams);
					if(pParams->nAlgorithm == _N_ALGORITHM_STOCHASTIC_II)	/* lookup table of probabilities is applied before consWeight */
					{
						/* replace with value from lookup table */
						lookupPos = (int)(pLandscape->asUndev[pLandscape->undevSites].logitVal * (pParams->nProbLookup - 1));
						pLandscape->asUndev[pLandscape->undevSites].logitVal = pParams->adProbLookup[lookupPos];
//						fprintf(stdout, "%f %d %f\n", pLandscape->asUndev[pLandscape->undevSites].logitVal, lookupPos, pParams->adProbLookup[lookupPos]);
					}
					//multiply consweight
					pLandscape->asUndev[pLandscape->undevSites].logitVal *= pThis->consWeight;
					pLandscape->asUndev[pLandscape->undevSites].bUntouched = pThis->bUntouched;	/* need to store this to put correct elements near top of list */
					if(pLandscape->asUndev[pLandscape->undevSites].logitVal > 0.0)
					{
						/* only add one more to the list if have a positive probability */
						pLandscape->undevSites++;
					}
				}
			}
		}
	}
	/* downweight the devPressure if necessary (do not do in first step) */
	/* doing it here means that last time step values have full weight */
	if(pParams->nAlgorithm == _N_ALGORITHM_STOCHASTIC_I)
	{
		if(pParams->dDevPersistence < 1.0)
		{
			fprintf(stdout, "\t\tdownweighting development pressure\n");

			for(i=0;i<pLandscape->totalCells;i++)
			{
				pThis = &(pLandscape->asCells[i]);
				if(pThis->nCellType ==_CELL_VALID)
				{
					if(pThis->bUndeveloped)	/* only need to bother downweighting on cells that can still convert */
					{
						pThis->devPressure = (int) ((double)pThis->devPressure * pParams->dDevPersistence);
					}
				}
			}
		}
	}
	/* sort */
	if(pParams->sortProbs)	/* can only be zero for algorithm=stochastic_ii */
	{
		fprintf(stdout, "\t\tsorting %d unconserved undeveloped sites\n", pLandscape->undevSites);
		qsort(pLandscape->asUndev,pLandscape->undevSites,sizeof(t_Undev),undevCmpReverse);
	}
	else
	{
		fprintf(stdout, "\t\tskipping sort as choosing cells randomly\n");
	}
	//calculate cumulative probability // From Wenwu Tang
    double sum=pLandscape->asUndev[0].logitVal;
    for(i=1;i<pLandscape->undevSites;i++){
        pLandscape->asUndev[i].cumulProb=pLandscape->asUndev[i-1].cumulProb+pLandscape->asUndev[i].logitVal;
        sum=sum+pLandscape->asUndev[i].logitVal;
    }
    for(i=0;i<pLandscape->undevSites;i++){
        pLandscape->asUndev[i].cumulProb=pLandscape->asUndev[i].cumulProb/sum;
    }
}

/*
	depreciated method of creating patches based on spiraling around
*/
int fillValidNeighbourList(int nThisID, t_Landscape *pLandscape, t_Params *pParams, int *anToConvert, int nWantToConvert, int bAllowTouched, int bDeterministic)
{
	int nTried,nFound,x,y,upDown,stopAt,countMove,stepMove,nPos,nSkipped;

	anToConvert[0] = nThisID;
	nSkipped = 0;
	nFound = 1;
	nTried = 1;
	x = pLandscape->asCells[nThisID].thisX;
	y = pLandscape->asCells[nThisID].thisY;
	stopAt = 0;
	upDown = 0;
	stepMove = -1;
	while(nFound < nWantToConvert && ((double)nWantToConvert * pParams->giveUpRatio) > nTried)
	{
		countMove = 0;
		upDown = !upDown;
		if(upDown)
		{
			stopAt++;
			stepMove *= -1;
		}
		while(countMove < stopAt && nFound < nWantToConvert && ((double)nWantToConvert  * pParams->giveUpRatio) > nTried)
		{
			if(upDown)
			{
				x += stepMove;
			}
			else
			{
				y += stepMove;
			}
			//fprintf(stdout, "%d %d\n", x, y);
			nPos = posFromXY(x,y, pLandscape);
			if(nPos != _CELL_OUT_OF_RANGE)
			{
				if(pLandscape->asCells[nPos].nCellType == _CELL_VALID)
				{
					if(pLandscape->asCells[nPos].bUndeveloped)
					{
						if(bAllowTouched || pLandscape->asCells[nPos].bUntouched)
						{
							if(bDeterministic || (uniformRandom() < pParams->dProbWeight))
							{
								if(pLandscape->asCells[nPos].consWeight > 0.0)
								{
									anToConvert[nFound] = nPos;
									nFound++;
								}
							}
							else
							{
								pLandscape->asCells[nPos].bUntouched = 0;
								nSkipped++;
							}
						}
					}
				}
			}
			nTried++;
			countMove++;
		}
	}
//	fprintf(stdout, "\tskipped=%d\n", nSkipped);
	return nFound;
}

/*
	helper structures for neighbour grabbing algorithm
*/
typedef struct
{
	double				probAdd;
	int					cellID;
	int					newInList;
	double              distance;// distance to the center //Wenwu Tang
} t_candidateNeighbour;

typedef struct
{
	double					maxProb;
	int						nCandidates;
	int						nSpace;
	t_candidateNeighbour 	*aCandidates;
} t_neighbourList;

#define 	_N_NEIGHBOUR_LIST_BLOCK_SIZE	20

int addNeighbourIfPoss(int x, int y, t_Landscape *pLandscape, t_neighbourList *pNeighbours, t_Params *pParams)
{
	int 	i,thisPos,mustAdd,listChanged,lookupPos;
	double	probAdd;
	t_Cell  *pThis;

	listChanged=0;
	thisPos = posFromXY(x, y, pLandscape);
	if(thisPos != _CELL_OUT_OF_RANGE)
	{
		pThis = &(pLandscape->asCells[thisPos]);
		if(pThis->nCellType == _CELL_VALID)
		{
			pThis->bUntouched = 0;
			if(pThis->bUndeveloped)
			{
				if(pThis->consWeight>0.0)
				{
					/* need to add this cell... */

					/* ...either refresh its element in list if already there */
					mustAdd = 1;
					for(i=0;mustAdd && i<pNeighbours->nCandidates;i++)
					{
						if(pNeighbours->aCandidates[i].cellID == thisPos)
						{
							pNeighbours->aCandidates[i].newInList=1;
							mustAdd=0;
							listChanged=1;
						}
					}
					/* or add it on the end, allocating space if necessary */
					if(mustAdd)
					{
						if(pNeighbours->nCandidates == pNeighbours->nSpace)
						{
							pNeighbours->nSpace+=_N_NEIGHBOUR_LIST_BLOCK_SIZE;
							pNeighbours->aCandidates = (t_candidateNeighbour*)realloc(pNeighbours->aCandidates, pNeighbours->nSpace*sizeof(t_candidateNeighbour));
							if(!pNeighbours->aCandidates)
							{
								fprintf(stderr, "memory error in addNeighbourIfPoss()\n...exiting\n");
								exit(EXIT_FAILURE);
							}
						}
						pNeighbours->aCandidates[pNeighbours->nCandidates].cellID = thisPos;
						pNeighbours->aCandidates[pNeighbours->nCandidates].newInList = 1;

						/* note duplication of effort in here recalculating the probabilities, but didn't store them in an accessible way */
						/*
						probAdd = pParams->dIntercept;
						probAdd += pParams->dDevPressure * pThis->devPressure;
						probAdd += pParams->dEmployAttraction * pThis->employAttraction;
						probAdd += pParams->dInterchangeDistance * pThis->interchangeDistance;
						probAdd += pParams->dRoadDensity * pThis->roadDensity;
						probAdd = 1.0/(1.0 + exp(probAdd));
                        */
						probAdd=getDevProbability(pThis,pParams);
						/* replace with value from lookup table */
						lookupPos = (int)(probAdd * (pParams->nProbLookup - 1));
						probAdd = pParams->adProbLookup[lookupPos];
						probAdd *= pThis->consWeight;
						/* multiply by the patch factor to help control shapes of patches */
						probAdd *= pParams->patchFactor;// pF now is set to 1, so it won't affect the result.or this can be deleted//WT
						pNeighbours->aCandidates[pNeighbours->nCandidates].probAdd = probAdd;
						/* only actually add it if will ever transition */
						if(probAdd > 0.0)
						{
							pNeighbours->nCandidates++;
							if(probAdd > pNeighbours->maxProb)
							{
								pNeighbours->maxProb = probAdd;
							}
							listChanged=1;
						}
					}
				}
			}
		}
	}
	return listChanged;
}
//From Wenwu: sort according to the new flag and conversion probability
static int sortNeighbours(const void *pOne, const void *pTwo)
{
	t_candidateNeighbour	*pNOne = (t_candidateNeighbour *)pOne;
	t_candidateNeighbour	*pNTwo = (t_candidateNeighbour *)pTwo;
/*
	if(pNOne->newInList > pNTwo->newInList)
	{
		return -1;
	}
	if(pNTwo->newInList > pNOne->newInList)
	{
		return 1;
	}
*/
	float p=rand()/(float)RAND_MAX;
	/*
	if(p<0.0){
	if(pNOne->distance> pNTwo->distance)
	{
		return 1;
	}
	if(pNTwo->distance > pNOne->distance)
	{
		return -1;
	}
    }
    */
    float alpha=0.5; // assign this from parameter file
    alpha=(sParams.patchMean)-(sParams.patchRange)*0.5;
    alpha+=p*sParams.patchRange;
	//alpha=0;	
	float p1,p2;
	p1=pNOne->probAdd/pow(pNOne->distance,alpha);
	p2=pNTwo->probAdd/pow(pNTwo->distance,alpha);

	//p1=alpha*pNOne->probAdd+(1-alpha)*pow(pNOne->distance,-2);
	//p2=alpha*pNTwo->probAdd+(1-alpha)*pow(pNTwo->distance,-2);

	if(p1 > p2)
	{
		return -1;
	}
	if(p2> p1)
	{
		return 1;
	}
	return 0;
}
//randomly shuffle the 4 neighbors
void Shuffle4(int*mat){
    int proba[4],flag[4];
    int size=4,i;
    for(i=0;i<size;i++){
        proba[i]=rand()%1000;
        flag[i]=1;
    }

    int numChecked=0;
    int max; int index=0;
    while(numChecked!=4){
        max=-9999;
        for(i=0;i<size;i++){
            if(flag[i]!=0){
                    if(proba[i]>max){
                        max=proba[i];
                        index=i;
                    }
            }
        }
        flag[index]=0;
        mat[numChecked]=index;
        numChecked++;
    }
}
void findNeighbours(int nThisID, t_Landscape *pLandscape, t_neighbourList *pNeighbours, t_Params *pParams)
{
	t_Cell	*pThis;
	int		listChanged=0;
    int     idmat[4]; idmat[0]=0;idmat[1]=1;idmat[2]=2;idmat[3]=3;
    int     i;
	/* try to add the neighbours of this one */
	pThis = &(pLandscape->asCells[nThisID]);
    //note: if sorted, then shuffle is no need
    /*
    Shuffle4(&idmat[0]);
    for(i=0;i<4;i++){
        if(idmat[i]==0)
            listChanged+=addNeighbourIfPoss(pThis->thisX-1,pThis->thisY,pLandscape,pNeighbours,pParams);
        else if(idmat[i]==1)
            listChanged+=addNeighbourIfPoss(pThis->thisX+1,pThis->thisY,pLandscape,pNeighbours,pParams);
        else if(idmat[i]==2)
            listChanged+=addNeighbourIfPoss(pThis->thisX,pThis->thisY-1,pLandscape,pNeighbours,pParams);
        else if (idmat[i]==3)
            listChanged+=addNeighbourIfPoss(pThis->thisX,pThis->thisY+1,pLandscape,pNeighbours,pParams);
    }
    */

    listChanged+=addNeighbourIfPoss(pThis->thisX-1,pThis->thisY,pLandscape,pNeighbours,pParams);//left
    listChanged+=addNeighbourIfPoss(pThis->thisX+1,pThis->thisY,pLandscape,pNeighbours,pParams);//right
    listChanged+=addNeighbourIfPoss(pThis->thisX,pThis->thisY-1,pLandscape,pNeighbours,pParams);//down
    listChanged+=addNeighbourIfPoss(pThis->thisX,pThis->thisY+1,pLandscape,pNeighbours,pParams);//up 
    if(sParams.numNeighbors==8){
        listChanged+=addNeighbourIfPoss(pThis->thisX-1,pThis->thisY-1,pLandscape,pNeighbours,pParams);
        listChanged+=addNeighbourIfPoss(pThis->thisX-1,pThis->thisY+1,pLandscape,pNeighbours,pParams);
        listChanged+=addNeighbourIfPoss(pThis->thisX+1,pThis->thisY-1,pLandscape,pNeighbours,pParams);
        listChanged+=addNeighbourIfPoss(pThis->thisX+1,pThis->thisY+1,pLandscape,pNeighbours,pParams);
    }
	/* if anything has been altered then resort */

	if(listChanged)
	{
		//qsort(pNeighbours->aCandidates,pNeighbours->nCandidates,sizeof(t_candidateNeighbour),sortNeighbours);
	}
}
double getDistance(int id1,int id2,t_Landscape *pLandscape){
    double x1,x2,y1,y2, result=0;
    x1=pLandscape->asCells[id1].thisX;
    x2=pLandscape->asCells[id2].thisX;
    y1=pLandscape->asCells[id1].thisY;
    y2=pLandscape->asCells[id2].thisY;
    result=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return result;
}
double getDistance1(double x1,double y1,double x2,double y2){
    double result=0;
    result=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    return result;
}
int newPatchFinder(int nThisID, t_Landscape *pLandscape, t_Params *pParams, int *anToConvert, int nWantToConvert)
{
	int 			bTrying,i,nFound,j;
	t_neighbourList	sNeighbours;

	memset(&sNeighbours,0,sizeof(t_neighbourList));
	anToConvert[0] = nThisID;
	pLandscape->asCells[nThisID].bUntouched = 0;
	pLandscape->asCells[nThisID].bUndeveloped = 0;
	nFound=1;
	if(pLandscape->asCells[nThisID].index_region==2)
        int stop=1;
	findNeighbours(nThisID,pLandscape,&sNeighbours,pParams);
	for(i=0;i<sNeighbours.nCandidates;i++){
        sNeighbours.aCandidates[i].distance=getDistance(nThisID,sNeighbours.aCandidates[i].cellID,pLandscape);
	}
	//print2ASC(pLandscape,"./test.txt");
	while(nFound < nWantToConvert && sNeighbours.nCandidates > 0)
	{
		i=0;
		bTrying=1;
		while(bTrying)
		{
			if(uniformRandom() < sNeighbours.aCandidates[i].probAdd)
			{
				anToConvert[nFound] = sNeighbours.aCandidates[i].cellID;
				/* flag that it is about to develop */
				pLandscape->asCells[anToConvert[nFound]].bUndeveloped = 0;
				/* remove this one from the list by copying down everything above it */
				for(j=i+1;j<sNeighbours.nCandidates;j++)
				{
					sNeighbours.aCandidates[j-1].cellID = sNeighbours.aCandidates[j].cellID;
					sNeighbours.aCandidates[j-1].newInList = sNeighbours.aCandidates[j].newInList;
					sNeighbours.aCandidates[j-1].probAdd = sNeighbours.aCandidates[j].probAdd;
				}
				/* reduce the size of the list */
				sNeighbours.nCandidates--;
				sNeighbours.maxProb = 0.0;
				for(j=0;j<sNeighbours.nCandidates;j++)
				{
					if(sNeighbours.aCandidates[j].probAdd > sNeighbours.maxProb)
					{
						sNeighbours.maxProb = sNeighbours.aCandidates[j].probAdd;
					}
				}
				/* add its neighbours */
				findNeighbours(anToConvert[nFound],pLandscape,&sNeighbours,pParams);
				nFound++;
				bTrying=0;
			}
			else
			{
				sNeighbours.aCandidates[i].newInList=0; 
			}
			if(bTrying)
			{
				i++;
				if(i==sNeighbours.nCandidates)
				{
					i=0;
				}
			}
			else
			{
				/* always resort the list if have just added, to let new elements bubble to the top */
				if(sNeighbours.nCandidates>0)
				{
//					int z;
                    for(i=0;i<sNeighbours.nCandidates;i++){
                        sNeighbours.aCandidates[i].distance=getDistance(nThisID,sNeighbours.aCandidates[i].cellID,pLandscape);
                    }

					qsort(sNeighbours.aCandidates,sNeighbours.nCandidates,sizeof(t_candidateNeighbour),sortNeighbours);
#if 0
					for(z=0;z<sNeighbours.nCandidates;z++)
					{
						fprintf(stdout, "%d %d %f\n", z, sNeighbours.aCandidates[z].newInList, sNeighbours.aCandidates[z].probAdd);
					}
					fprintf(stdout, "***\n");
#endif
				}
			}
		}
	}
	if(sNeighbours.nSpace)	/* free any allocated memory */
	{
		free(sNeighbours.aCandidates);
	}
//	fprintf(stdout, "looking for %d found %d\n", nWantToConvert,nFound);
	return nFound;
}

int convertCells(t_Landscape *pLandscape, t_Params *pParams, int nThisID, int nStep, int bAllowTouched, int bDeterministic)
{
	int			i,x,y,nCell,nDone,nToConvert,nWantToConvert;
	int			*anToConvert;
	t_Cell		*pThis,*pThat;

	nDone = 0;
	/* in here need to build list of near neigbours to loop over */
	nWantToConvert = pLandscape->aParcelSizes[(int)(uniformRandom()*pLandscape->parcelSizes)];
	anToConvert = (int*)malloc(sizeof(int) * nWantToConvert);
	if(anToConvert)
	{
		/* in here goes code to fill up list of neighbours */
		if(pParams->nAlgorithm == _N_ALGORITHM_STOCHASTIC_II)
		{
        nToConvert = newPatchFinder(nThisID,pLandscape,pParams,anToConvert, nWantToConvert);
		}
		else
		{
			nToConvert = fillValidNeighbourList(nThisID, pLandscape, pParams, anToConvert, nWantToConvert, bAllowTouched, bDeterministic);
		}
		if(nToConvert > 0)	/* will actually always be the case */
		{
//			fprintf(stdout, "wanted %d, found %d\n", nWantToConvert, nToConvert);
			for(i=0;i<nToConvert;i++) 
			{
				/* convert, updating dev press on neighbours */
				pThis = &(pLandscape->asCells[anToConvert[i]]);
				pThis->bUndeveloped = 0;
				pThis->bUntouched = 0;
				pThis->tDeveloped = nStep+1;
				nDone++;
				// fprintf(stdout, "%d total=%f dev=%f employment=%f interchange=%f roads=%f\n", pLandscape->asUndev[i].cellID, pLandscape->asUndev[i].logitVal, pParams->dDevPressure * pThis->devPressure, pParams->dEmployAttraction * pThis->employAttraction, pParams->dInterchangeDistance * pThis->interchangeDistance, pParams->dRoadDensity * pThis->roadDensity);
				float dist=0;
				float force=0;
				for(x = pThis->thisX - pParams->nDevNeighbourhood;x <= pThis->thisX + pParams->nDevNeighbourhood; x++)
				{
					for(y = pThis->thisY - pParams->nDevNeighbourhood;y <= pThis->thisY + pParams->nDevNeighbourhood; y++)
					{
						nCell = posFromXY(x,y,pLandscape);
						if(nCell != _CELL_OUT_OF_RANGE)
						{
							pThat = &(pLandscape->asCells[nCell]);
							dist=getDistance1(pThis->thisX,pThis->thisY,x,y);
							if(dist>pParams->nDevNeighbourhood) continue;
                            force=0;
							if(pParams->devPressureApproach==1)
                                force=1; // development pressure increases by 1
                            else if(pParams->devPressureApproach==2){
                                if(dist>0) force=pParams->scalingFactor/pow(dist,pParams->alpha);
                            }
                            else{
                                force=pParams->scalingFactor*exp(-2*dist/pParams->alpha);
                            }
                            pThat->devPressure=pThat->devPressure+force;
						}
					}
				}
			}
		}
		free(anToConvert);
	}
	return nDone;
}

/*
	main routine to actually run the model
*/
void updateMap(t_Landscape *pLandscape, t_Params *pParams)
{
	FILE		*fIn;
	char		*szBuff;
	char		szStepLabel[_N_MAX_FILENAME_LEN];
	int			i,nToConvert,nStep,nDone,bAllowTouched,nRandTries,nExtra;
	double		dProb;
	t_Cell		*pThis;

	nExtra=0;
	fprintf(stdout, "entered updateMap()\n");
	fIn = fopen(pParams->controlFile, "rb");
	if(fIn)
	{
		szBuff = (char*)malloc(_N_MAX_DYNAMIC_BUFF_LEN*sizeof(char));
		if(szBuff)
		{
			nStep = 1;	/* start counter at 1 so can distinguish between cells which transition on first step and those which already developed */
			while(fgets(szBuff, _N_MAX_DYNAMIC_BUFF_LEN, fIn))
			{
				if(parseControlLine(szBuff,szStepLabel,&nToConvert))
				{
					fprintf(stdout, "\tdoing step %s...controlFile requests conversion of %d cells\n", szStepLabel, nToConvert);
					if(nExtra > 0)
					{
						if(nToConvert - nExtra > 0)
						{
							nToConvert -= nExtra;
							nExtra = 0;
						}
						else
						{
							nExtra -= nToConvert;
							nToConvert = 0;
						}
					}
					fprintf(stdout, "\t\tafter accounting for extra cells, attempt %d cells\n", nToConvert);
					/* if have cells to convert this step */
					if(nToConvert > 0)
					{
						findAndSortProbs(pLandscape, pParams, nToConvert);
						/* if not enough cells to convert then alter number required */
						if(nToConvert > pLandscape->undevSites)
						{
							fprintf(stdout, "\t\tnot enough undeveloped sites...converting all\n");
							nToConvert = pLandscape->undevSites;
						}
						/* update either in deterministic or stochastic fashion */
						fprintf(stdout, "\t\tupdating map\n");
						switch(pParams->nAlgorithm)
						{
						case _N_ALGORITHM_DETERMINISTIC:	/* deterministic */
							nDone = 0;
							for(i=0;i<nToConvert&&nDone<nToConvert;i++)
							{
								nDone+=convertCells(pLandscape, pParams, pLandscape->asUndev[i].cellID, nStep, 1, 1);
							}
							break;
						case _N_ALGORITHM_STOCHASTIC_I: /* stochastic */
							nDone = 0;
							i = 0;
							bAllowTouched = 0;
							while(nDone < nToConvert)	/* loop until done enough cells...might need multiple passes */
							{
								if(i == pLandscape->undevSites)
								{
									i = 0;					/* if at the end of the grid, just loop around again until done */
									bAllowTouched = 1;		/* allow previously considered cells if you have to */
								}
								pThis = &(pLandscape->asCells[pLandscape->asUndev[i].cellID]);
								if(pThis->bUndeveloped)		/* need to check is still undeveloped */
								{
									if(bAllowTouched || pThis->bUntouched)
									{
										/* Doug's "back to front" logit */
//										dProb = 1.0/(1.0 + exp(pLandscape->asUndev[i].logitVal));
										dProb = pLandscape->asUndev[i].logitVal;
										/* if starting a patch off here */
										if(uniformRandom() < pParams->dProbWeight * dProb)
										{
											nDone+=convertCells(pLandscape, pParams, pLandscape->asUndev[i].cellID, nStep, bAllowTouched, 0);
										}
									}
									pThis->bUntouched = 0;
								}
								i++;
							}
							break;
						case _N_ALGORITHM_STOCHASTIC_II: /* stochastic */
							nDone = 0;
							i = 0;
							nRandTries = 0;
							bAllowTouched = 0;
							while(nDone < nToConvert)	/* loop until done enough cells...might need multiple passes */
							{
								if(i == pLandscape->undevSites)
								{
									i = 0;					/* if at the end of the grid, just loop around again until done */
									bAllowTouched = 1;		/* allow previously considered cells if you have to */
								}
								/* if tried too many randomly in this step, give up on idea of only letting untouched cells convert */
								if(nRandTries > _MAX_RAND_FIND_SEED_FACTOR * nToConvert )
								{
									bAllowTouched = 1;
								}
								if(pParams->sortProbs)
								{
									/* if sorted then choose the top cell and do nothing */
								}
								else
								{
									/* otherwise give a random undeveloped cell a go */
                                   if(sParams.seedSearch==1)
                                        i = (int)(uniformRandom()*pLandscape->undevSites);
									//pick one according to their probability
                                   else
                                        i=getUnDevIndex(pLandscape);

								}
								pThis = &(pLandscape->asCells[pLandscape->asUndev[i].cellID]);
								//print2ASC(pLandscape,"./test.txt");
								if(pThis->bUndeveloped)		/* need to check is still undeveloped */
								{
									if(bAllowTouched || pThis->bUntouched)
									{
										/* Doug's "back to front" logit */
//										dProb = 1.0/(1.0 + exp(pLandscape->asUndev[i].logitVal));
										dProb = pLandscape->asUndev[i].logitVal;
										if(uniformRandom() < dProb)
										{
											nDone+=convertCells(pLandscape, pParams, pLandscape->asUndev[i].cellID, nStep, bAllowTouched, 0);
										}
									}
									pThis->bUntouched = 0;
								}
								if(pParams->sortProbs)
								{
									i++;
								}
								else
								{
									nRandTries++;
								}
							}
							break;
						default:
							fprintf(stderr, "Unknown algorithm...exiting\n");
							break;
						}
						fprintf(stdout, "\t\tconverted %d sites\n", nDone);
						nExtra += (nDone - nToConvert);
						fprintf(stdout, "\t\t%d extra sites knocked off next timestep\n", nExtra);
					}
				}
				nStep++;	/* next time step */
			}
			/* dump results to a file at the end of the run */
			dumpDevRaster(pLandscape, pParams, pParams->dumpFile);
			free(szBuff);
		}
		fclose(fIn);
	}
}

void readDevDemand(t_Params *pParams){
    ifstream f1;//this one may be deleted
/* 
   f1.open(pParams->controlFile);
    int counter=0;
    int val1,val2;
    while(!f1.eof()){
        f1>>val1>>val2;
        pParams->devDemand[counter]=val2;
        counter++;
    }
    pParams->nSteps=counter;
    f1.close();
*/  
  //read land demand for each region
    int i;
    int val1,val2,counter;
	char str[200];
    f1.open(pParams->controlFileAll);
    f1>>str>>val1;
    pParams->nSteps=val1;
    f1.getline(str,200);
    f1.getline(str,200); //skip the header
    counter=0; int ii;
    for(ii=0;ii<pParams->nSteps;ii++){
        f1>>val1;
        for(i=0;i<pParams->num_Regions;i++){
            f1>>val1;
            pParams->devDemands[i][ii]=val1;
        }

    }
    f1.close();

}
void initializeUnDev(t_Landscape *pLandscape, t_Params *pParams){
    int i;
    pLandscape->asUndevs=new t_Undev*[pParams->num_Regions];
    for(i=0;i<pParams->num_Regions;i++){
        pLandscape->asUndevs[i]=new t_Undev[MAX_UNDEV_SIZE];
        pLandscape->num_undevSites[i]=0;
    }
}
void finalizeUnDev(t_Landscape *pLandscape, t_Params *pParams){
    int i;

    for(i=0;i<pParams->num_Regions;i++){
        delete [] pLandscape->asUndevs[i];
    }
    delete [] pLandscape->asUndevs;
}
/*
	main routine to run models on multiple regions //WTang
*/
void updateMapAll(t_Landscape *pLandscape, t_Params *pParams){
    //readDevDemand(pParams);
    int i,j;
    //initialize undev arrays
    initializeUnDev(pLandscape, pParams);

    //iterate each step (year)
    for(i=0;i<pParams->nSteps;i++){
        cout<<i<<"\t"<<pParams->nSteps<<endl;
        if(i==10)
            int stop=1;
        findAndSortProbsAll(pLandscape, pParams,i);//for each sub-region, find and update conversion probability (conservation weight applied)
        for(j=0;j<pParams->num_Regions;j++)
	//j=1;
            updateMap1(pLandscape, pParams, i,j);
    }
    dumpDevRaster(pLandscape, pParams, pParams->dumpFile);
    finalizeUnDev(pLandscape, pParams);
}

void updateMap1(t_Landscape *pLandscape, t_Params *pParams, int step, int regionID)
{
	int			i,nToConvert,nStep,nDone,bAllowTouched,nRandTries,nExtra;
	double		dProb;
	t_Cell		*pThis;
	nExtra=0;

    nStep=step;

	//fprintf(stdout, "entered updateMap()\n");
    nToConvert=pParams->devDemands[regionID][step];
    //fprintf(stdout, "\tdoing step %s...controlFile requests conversion of %d cells\n", szStepLabel, nToConvert);
    if(nExtra > 0){
        if(nToConvert - nExtra > 0){
            nToConvert -= nExtra;
            nExtra = 0;
        }
        else{
            nExtra -= nToConvert;
			nToConvert = 0;
        }
    }
    fprintf(stdout, "\t\tafter accounting for extra cells, attempt %d cells\n", nToConvert);
    /* if have cells to convert this step */
    if(nToConvert > 0){
        //findAndSortProbs(pLandscape, pParams, nToConvert);
        /* if not enough cells to convert then alter number required */
		if(nToConvert > pLandscape->num_undevSites[regionID]){
            fprintf(stdout, "\t\tnot enough undeveloped sites...converting all\n");
			nToConvert = pLandscape->num_undevSites[regionID];
        }
		/* update either in deterministic or stochastic fashion */
		//fprintf(stdout, "\t\tupdating map\n");
		switch(pParams->nAlgorithm)
		{
            case _N_ALGORITHM_DETERMINISTIC:	/* deterministic */
                nDone = 0;
				for(i=0;i<nToConvert&&nDone<nToConvert;i++)
				{
					nDone+=convertCells(pLandscape, pParams, pLandscape->asUndev[i].cellID, nStep, 1, 1);
                }
				break;
            case _N_ALGORITHM_STOCHASTIC_I: /* stochastic */
				nDone = 0;
				i = 0;
				bAllowTouched = 0;
				while(nDone < nToConvert)	/* loop until done enough cells...might need multiple passes */
				{
                    if(i == pLandscape->undevSites)
					{
                        i = 0;					/* if at the end of the grid, just loop around again until done */
						bAllowTouched = 1;		/* allow previously considered cells if you have to */
                    }
					pThis = &(pLandscape->asCells[pLandscape->asUndev[i].cellID]);
					if(pThis->bUndeveloped)		/* need to check is still undeveloped */
					{
                        if(bAllowTouched || pThis->bUntouched)
						{
                            /* Doug's "back to front" logit */
    						//dProb = 1.0/(1.0 + exp(pLandscape->asUndev[i].logitVal));
                            dProb = pLandscape->asUndev[i].logitVal;
                            /* if starting a patch off here */
							if(uniformRandom() < pParams->dProbWeight * dProb)
							{
                                nDone+=convertCells(pLandscape, pParams, pLandscape->asUndev[i].cellID, nStep, bAllowTouched, 0);
                            }
                        }
						pThis->bUntouched = 0;
                    }
                    i++;
				}
            break;
            case _N_ALGORITHM_STOCHASTIC_II: /* stochastic */
                nDone = 0;
				i = 0;
				nRandTries = 0;
				bAllowTouched = 0;
				while(nDone < nToConvert)	/* loop until done enough cells...might need multiple passes */
				{
                    if(i == pLandscape->num_undevSites[regionID])
					{
                        i = 0;					/* if at the end of the grid, just loop around again until done */
						bAllowTouched = 1;		/* allow previously considered cells if you have to */
                    }
					/* if tried too many randomly in this step, give up on idea of only letting untouched cells convert */
					if(nRandTries > _MAX_RAND_FIND_SEED_FACTOR * nToConvert )
					{
                        bAllowTouched = 1;
                    }
					if(pParams->sortProbs)
                    {
                        /* if sorted then choose the top cell and do nothing */
                    }
					else
					{
                        /* otherwise give a random undeveloped cell a go */
                        if(sParams.seedSearch==1)
                            i = (int)(uniformRandom()*pLandscape->num_undevSites[regionID]);
							//pick one according to their probability
                        else
                            i=getUnDevIndex1(pLandscape,regionID);
                    }
					pThis = &(pLandscape->asCells[pLandscape->asUndevs[regionID][i].cellID]);
					//print2ASC(pLandscape,"./test.txt");
					if(pThis->bUndeveloped)		/* need to check is still undeveloped */
					{
                        if(bAllowTouched || pThis->bUntouched)
						{
                            /* Doug's "back to front" logit */
//							dProb = 1.0/(1.0 + exp(pLandscape->asUndev[i].logitVal));
							dProb = pLandscape->asUndevs[regionID][i].logitVal;
							if(uniformRandom() < dProb)
							{
                                nDone+=convertCells(pLandscape, pParams, pLandscape->asUndevs[regionID][i].cellID, nStep, bAllowTouched, 0);
                            }
                        }
						pThis->bUntouched = 0;
                    }
					if(pParams->sortProbs)
					{
                        i++;
                    }
					else
					{
                        nRandTries++;
                    }
                }
                break;
				default:
                    fprintf(stderr, "Unknown algorithm...exiting\n");
					break;
        }
		fprintf(stdout, "\t\tconverted %d sites\n", nDone);
		nExtra += (nDone - nToConvert);
		fprintf(stdout, "\t\t%d extra sites knocked off next timestep\n", nExtra);
    }
}


/*
	check Doug's calculation of devPressure...no longer called
*/
void	testDevPressure(t_Landscape *pLandscape, t_Params *pParams)
{
	int x,y,xN,yN,cellID,nCell,testPressure;

	for(y=0;y<pLandscape->maxY;y++)
	{
		for(x=0;x<pLandscape->maxX;x++)
		{
			cellID = posFromXY(x, y, pLandscape);
			if(cellID != _CELL_OUT_OF_RANGE)	/* should never happen */
			{
				if(pLandscape->asCells[cellID].nCellType == _CELL_VALID)
				{
					if(pLandscape->asCells[cellID].bUndeveloped)
					{
						if(pLandscape->asCells[cellID].devPressure)
						{
							fprintf(stdout, "(x,y)=(%d,%d) cellType=%d bUndev=%d devPress=%d\n", x, y, pLandscape->asCells[cellID].nCellType, pLandscape->asCells[cellID].bUndeveloped, pLandscape->asCells[cellID].devPressure);
							testPressure = 0;
							for(xN = x - pParams->nDevNeighbourhood;xN <= x + pParams->nDevNeighbourhood; xN++)
							{
								for(yN = y - pParams->nDevNeighbourhood;yN <= y + pParams->nDevNeighbourhood; yN++)
								{
									nCell = posFromXY(xN,yN,pLandscape);
									if(nCell != _CELL_OUT_OF_RANGE)
									{
										if(pLandscape->asCells[nCell].nCellType == _CELL_VALID)
										{
											if(pLandscape->asCells[nCell].bUndeveloped == 0)
											{
												testPressure++;
											}
										}
									}
								}
							}
							fprintf(stdout, "\t%d\n", testPressure);
						}
					}
				}
			}
		}
	}
}

int readParcelSizes(t_Landscape *pLandscape, t_Params *pParams)
{
	FILE	*fIn;
	char	*szBuff;
	int		nMaxParcels;
       

	pLandscape->parcelSizes = 0;

	fprintf(stdout, "entered readParcelSizes()\n");
	fIn = fopen(pParams->parcelSizeFile, "rb");
	if(fIn)
	{
		szBuff = (char*)malloc(_N_MAX_DYNAMIC_BUFF_LEN*sizeof(char));
		if(szBuff)
		{
			/* just scan the file twice */
			nMaxParcels = 0;
			while(fgets(szBuff, _N_MAX_DYNAMIC_BUFF_LEN,fIn))
			{
				nMaxParcels++;
			}
			rewind(fIn);
			if(nMaxParcels)
			{
				pLandscape->aParcelSizes= (int*)malloc(sizeof(int) * nMaxParcels);
				if(pLandscape->aParcelSizes)
				{
					while(fgets(szBuff, _N_MAX_DYNAMIC_BUFF_LEN,fIn))
					{
						pLandscape->aParcelSizes[pLandscape->parcelSizes] = atoi(szBuff)*pParams->discountFactor;
						if(pLandscape->aParcelSizes[pLandscape->parcelSizes])
						{
							pLandscape->parcelSizes++;
						}
					}
				}
			}
			free(szBuff);
		}
		fclose(fIn);
	}
	return pLandscape->parcelSizes;
}

/*
	main code.

	no command line args, but need argv for getCfgFileName()
*/
int main(int argc, char **argv)
{

    struct
    {
        struct Option
                *xSize, *ySize,
                *controlFile, *employAttractionFile, *interchangeDistanceFile,
                *roadDensityFile, *undevelopedFile, *devPressureFile, *consWeightFile,
                *addVariableFiles, *nDevNeighbourhood, *dumpFile, *algorithm, *dProbWeight,
                *dDevPersistence, *parcelSizeFile, *discountFactor, *giveUpRatio;

    } opt;
    struct
    {
        struct Flag *a;
    } flag;

    G_gisinit(argv[0]);

    struct GModule *module = G_define_module();
    G_add_keyword(_("raster"));
    G_add_keyword(_("patch growing"));
    G_add_keyword(_("urban"));
    G_add_keyword(_("landscape"));
    G_add_keyword(_("modeling"));
    module->label = _("FUTURES");
    module->description = _("...");

    opt.xSize = G_define_option();
    opt.xSize->key = "xsize";
    opt.xSize->type = TYPE_INTEGER;
    opt.xSize->required = YES;
    opt.xSize->description = _("Size of the grids");

    opt.ySize = G_define_option();
    opt.ySize->key = "ysize";
    opt.ySize->type = TYPE_INTEGER;
    opt.ySize->required = YES;
    opt.ySize->description = _("Size of the grids");

    opt.controlFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.controlFile->key = "control_file";
    opt.controlFile->type = TYPE_INTEGER;
    opt.controlFile->required = YES;
    opt.controlFile->description = _("File containing information on how many cells to transition and when");

    opt.employAttractionFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.employAttractionFile->key = "employ_attraction";
    opt.employAttractionFile->required = YES;
    opt.employAttractionFile->description = _("Files containing the information to read in");

    opt.interchangeDistanceFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.interchangeDistanceFile->key = "interchange_distance";
    opt.interchangeDistanceFile->required = YES;
    opt.interchangeDistanceFile->description = _("Files containing the information to read in");

    opt.roadDensityFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.roadDensityFile->key = "road_density";
    opt.roadDensityFile->required = YES;
    opt.roadDensityFile->description = _("Files containing the information to read in");

    opt.undevelopedFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.undevelopedFile->key = "undeveloped";
    opt.undevelopedFile->required = YES;
    opt.undevelopedFile->description = _("Files containing the information to read in");

    opt.devPressureFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.devPressureFile->key = "development_pressure";
    opt.devPressureFile->required = YES;
    opt.devPressureFile->description = _("Files containing the information to read in");

    opt.consWeightFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.consWeightFile->key = "cons_weight";
    opt.consWeightFile->required = YES;
    opt.consWeightFile->description = _("Files containing the information to read in");

    opt.addVariableFiles = G_define_standard_option(G_OPT_F_INPUT);
    opt.addVariableFiles->key = "additional_variable_files";
    opt.addVariableFiles->required = YES;
    opt.addVariableFiles->multiple = YES;
    opt.addVariableFiles->description = _("additional variables");

    opt.nDevNeighbourhood = G_define_option();
    opt.nDevNeighbourhood->key = "n_dev_neighbourhood";
    opt.nDevNeighbourhood->type = TYPE_INTEGER;
    opt.nDevNeighbourhood->required = YES;
    opt.nDevNeighbourhood->description = _("Size of square used to recalculate development pressure");

    opt.dumpFile = G_define_standard_option(G_OPT_F_OUTPUT);
    opt.dumpFile->key = "dump_file";
    opt.dumpFile->required = YES;
    opt.dumpFile->description = _("Used in writing rasters");

    opt.algorithm = G_define_option();
    opt.algorithm->key = "algorithm";
    opt.algorithm->type = TYPE_STRING;
    opt.algorithm->required = YES;
    opt.algorithm->options = "deterministic,stochastic1,stochastic2";
    opt.algorithm->description = _("Parameters controlling the algorithm to use");

    opt.dProbWeight = G_define_option();
    opt.dProbWeight->key = "prob_weight";
    opt.dProbWeight->type = TYPE_DOUBLE;
    opt.dProbWeight->required = YES;
    opt.dProbWeight->description = _("Parameters controlling the algorithm to use");

    opt.dDevPersistence = G_define_option();
    opt.dDevPersistence->key = "dev_neighbourhood";
    opt.dDevPersistence->type = TYPE_DOUBLE;
    opt.dDevPersistence->required = YES;
    opt.dDevPersistence->description = _("Parameters controlling the algorithm to use");

    opt.parcelSizeFile = G_define_standard_option(G_OPT_F_INPUT);
    opt.parcelSizeFile->key = "parcel_size_file";
    opt.parcelSizeFile->type = TYPE_DOUBLE;
    opt.parcelSizeFile->required = YES;
    opt.parcelSizeFile->description = _("File containing information on the parcel size to use");

    opt.discountFactor = G_define_option();
    opt.discountFactor->key = "discount_factor";
    opt.discountFactor->type = TYPE_DOUBLE;
    opt.discountFactor->required = YES;
    opt.discountFactor->description = _("discount factor of patch size");

    opt.giveUpRatio = G_define_option();
    opt.giveUpRatio->key = "give_up_ratio";
    opt.giveUpRatio->type = TYPE_DOUBLE;
    opt.giveUpRatio->required = YES;
    opt.giveUpRatio->description = _("Give up ratio");

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

	struct timeval ttime;
	gettimeofday(&ttime,NULL);

    // TODO: move this back to local variables
	//t_Params	sParams;
	//t_Landscape	sLandscape;

    seedRandom(ttime);

    /* blank everything out */
    memset(&sParams,0,sizeof(t_Params));

    /* set up parameters */
    sParams.xSize = atoi(opt.xSize->answer);
    sParams.ySize = atoi(opt.ySize->answer);
    sParams.controlFile = opt.controlFile->answer;
    sParams.employAttractionFile = opt.employAttractionFile->answer;
    sParams.interchangeDistanceFile = opt.interchangeDistanceFile->answer;
    sParams.roadDensityFile = opt.roadDensityFile->answer;
    sParams.undevelopedFile = opt.undevelopedFile->answer;
    sParams.devPressureFile = opt.devPressureFile->answer;
    sParams.consWeightFile = opt.consWeightFile->answer;
    sParams.numAddVariables = 0;
    char** answer = opt.addVariableFiles->answers;
    while (answer) {
        sParams.addVariableFile[sParams.numAddVariables] = *answer;
        sParams.numAddVariables += 1;
        // TODO: dyn allocate file list
        ++answer;
    }
    sParams.nDevNeighbourhood = atof(opt.nDevNeighbourhood->answer);
    sParams.dumpFile = opt.dumpFile->answer;
    if (!strcmp(opt.algorithm->answer, "deterministic,stochastic1,stochastic2"))
        sParams.nAlgorithm = _N_ALGORITHM_DETERMINISTIC;
    else if (!strcmp(opt.algorithm->answer, "stochastic1"))
        sParams.nAlgorithm = _N_ALGORITHM_STOCHASTIC_I;
    else if (!strcmp(opt.algorithm->answer, "stochastic2"))
        sParams.nAlgorithm = _N_ALGORITHM_STOCHASTIC_II;

    sParams.dProbWeight = atof(opt.dProbWeight->answer);
    sParams.dDevPersistence = atof(opt.dDevPersistence->answer);

    sParams.parcelSizeFile = opt.parcelSizeFile->answer;

    sParams.discountFactor = atof(opt.discountFactor->answer);
    sParams.giveUpRatio = atof(opt.giveUpRatio->answer);

    // TODO: always the same?
    sParams.sortProbs = 1;

    // TODO: _N_ALGORITHM_STOCHASTIC_II

			readDevDemand(&sParams);
			/* allocate memory */
			if(buildLandscape(&sLandscape, &sParams))
			{
				/* read data */
				if(readData(&sLandscape, &sParams))
				{
				    readData4AdditionalVariables(&sLandscape, &sParams);
					readIndexData(&sLandscape, &sParams);
					if(readParcelSizes(&sLandscape, &sParams))
					{
						//testDevPressure(&sLandscape, &sParams);
						/* do calculation and dump result */
						//one region, enable the following
						//updateMap(&sLandscape, &sParams);

						//multiple regions
						updateMapAll(&sLandscape, &sParams);
					}
					else
					{
						fprintf(stderr, "error in readParcelSizes()\n");
					}
				}
				else
				{
					fprintf(stderr, "error in readData()\n");
				}
				/* could put in routines to free memory, but OS will garbage collect anyway */
			}
			else
			{
				fprintf(stderr, "error in buildLandscape()\n");
			}
	//}//disable this when debuging

	//else
	//{
	//	fprintf(stderr, "error in getCfgFileName()...expecting file %s\n",szCfgFile);
	//}


	return EXIT_SUCCESS;
}

int getUnDevIndex(t_Landscape *pLandscape){
        float p=rand()/(double)RAND_MAX;
        int i;
        for(i=0;i<pLandscape->undevSites;i++){
            if(p<pLandscape->asUndev[i].cumulProb){
                return i;
            }
        }
}
int getUnDevIndex1(t_Landscape *pLandscape,int regionID){
        float p=rand()/(double)RAND_MAX;
        int i;
        for(i=0;i<pLandscape->num_undevSites[regionID];i++){
            if(p<pLandscape->asUndevs[regionID][i].cumulProb){
            	//cout<< "i "<< i<<" R "<<p<<", l "<<pLandscape->asUndevs[regionID][i].cumulProb<<endl;
                return i;
            }
        }
}

void print2ASC(t_Landscape*pLandscape,char* fn){
    ofstream f;
    f.open(fn);
    int i;

    for(i=0;i<pLandscape->maxX*pLandscape->maxY;i++){
        f<<pLandscape->asCells[i].bUndeveloped<<"\t";
        if(i%pLandscape->maxX==pLandscape->maxX-1)
            f<<endl;
    }
    f.close();
}

void findAndSortProbsAll(t_Landscape *pLandscape, t_Params *pParams,int step)
{
	int			i,lookupPos;
	t_Cell		*pThis;
    int         id;
	/* update calcs */
	fprintf(stdout, "\t\trecalculating probabilities\n");
    for(i=0;i<pParams->num_Regions;i++){
        pLandscape->num_undevSites[i]=0;
    }
	float val=0.0;
	for(i=0;i<pLandscape->totalCells;i++)
	{
		pThis = &(pLandscape->asCells[i]);
		if(pThis->nCellType ==_CELL_VALID)
		{
			if(pThis->bUndeveloped)
			{
				if(pThis->consWeight > 0.0)
				{
                    id=pThis->index_region-1;
                    //if(id>10)
                      //  int stop=1;
                    if(pThis->index_region==-9999) continue;
					/* note that are no longer just storing the logit value, but instead the probability (allows consWeight to affect sort order) */
					pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].cellID = i;
				val=getDevProbability(pThis,pParams);	
				pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].logitVal=val;
					pThis->devProba=val;
					if(pParams->nAlgorithm == _N_ALGORITHM_STOCHASTIC_II)	/* lookup table of probabilities is applied before consWeight */
					{
						/* replace with value from lookup table */
						lookupPos = (int)(pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].logitVal * (pParams->nProbLookup - 1));
						pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].logitVal = pParams->adProbLookup[lookupPos];
//						fprintf(stdout, "%f %d %f\n", pLandscape->asUndev[pLandscape->undevSites].logitVal, lookupPos, pParams->adProbLookup[lookupPos]);
					}
					pThis->devProba=pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].logitVal;
					pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].logitVal *= pThis->consWeight;//discount by a conservation factor
					pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].bUntouched = pThis->bUntouched;	/* need to store this to put correct elements near top of list */
					if(pLandscape->asUndevs[id][pLandscape->num_undevSites[id]].logitVal > 0.0)
					{
						/* only add one more to the list if have a positive probability */
						pLandscape->num_undevSites[id]++;
					}
				}
			}
		}
	}
	/* downweight the devPressure if necessary (do not do in first step) */
	/* doing it here means that last time step values have full weight */ 
	if(pParams->nAlgorithm == _N_ALGORITHM_STOCHASTIC_I)
	{
		if(pParams->dDevPersistence < 1.0)
		{
			fprintf(stdout, "\t\tdownweighting development pressure\n");

			for(i=0;i<pLandscape->totalCells;i++)
			{
				pThis = &(pLandscape->asCells[i]);
				if(pThis->nCellType ==_CELL_VALID)
				{
					if(pThis->bUndeveloped)	/* only need to bother downweighting on cells that can still convert */
					{
						pThis->devPressure = (int) ((double)pThis->devPressure * pParams->dDevPersistence);
					}
				}
			}
		}
	}
	/* sort */
	if(pParams->sortProbs)	/* can only be zero for algorithm=stochastic_ii */
	{
		fprintf(stdout, "\t\tsorting %d unconserved undeveloped sites\n", pLandscape->undevSites);
		qsort(pLandscape->asUndev,pLandscape->undevSites,sizeof(t_Undev),undevCmpReverse);
	}
	else
	{
		fprintf(stdout, "\t\tskipping sort as choosing cells randomly\n");
	}
	//calculate cumulative probability // From Wenwu Tang
	int j;
    for(j=0;j<pParams->num_Regions;j++){
        double sum=pLandscape->asUndevs[j][0].logitVal;
        for(i=1;i<pLandscape->num_undevSites[j];i++){
            pLandscape->asUndevs[j][i].cumulProb=pLandscape->asUndevs[j][i-1].cumulProb+pLandscape->asUndevs[j][i].logitVal;
            sum=sum+pLandscape->asUndevs[j][i].logitVal;
        }
        for(i=0;i<pLandscape->num_undevSites[j];i++){
            pLandscape->asUndevs[j][i].cumulProb=pLandscape->asUndevs[j][i].cumulProb/sum;
        }
    }
    char fn[200],str[20];
    strcpy(fn,"./proba");
    sprintf(str,"%d",step+1);
    strcat(fn,str); strcat(fn,".asc");
    //export2ASC1(pLandscape,pParams,0, fn);
}

void export2ASC1(t_Landscape *pLandscape, t_Params *pParams, int regionID,char* fn){
	int i,j;
	ofstream f;
	f.open(fn);
//	if(headflag==1){

	f<<"ncols	4008"<<endl;
	f<<"nrows	4678"<<endl;
//	f<<"xllcorner	516805"<<endl;
//	f<<"yllcorner	4945476"<<endl;
//	f<<"cellsize	100"<<endl;
	f<<"xllcorner	447482.1782041"<<endl;
	f<<"yllcorner	117225.9427485"<<endl;
	f<<"cellsize	30"<<endl;
	f<<"NODATA_value	-9999"<<endl;

//}
	t_Cell		*pThis;
    float val;
    int nc=pParams->xSize;
	for(i=0;i<pLandscape->totalCells;i++){
	    pThis = &(pLandscape->asCells[i]);
        val=pThis->devProba;
		f<<int(val*1000)<<"\t";
        if(i%nc==nc-1) f<<endl;
	}
	f.close();
}
