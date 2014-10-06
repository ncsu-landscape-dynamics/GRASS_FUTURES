#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "cfgutils.h"

/* work out configuration file name and check whether it exists */
int	getCfgFileName(char *szProgName, char *szCfgFile)
{
	char	*pPtr;
	FILE 	*fp;

	szCfgFile[0] = '\0';
	{
		if((pPtr = strrchr(szProgName,C_DIR_DELIMITER))!=NULL)
		{
			strcpy(szCfgFile,pPtr+1);
		}
		else
		{
			strcpy(szCfgFile,szProgName);
		}
		if((pPtr = strstr(szCfgFile,".exe"))!=NULL)
		{
			*pPtr = '\0';
		}
		strcat(szCfgFile,".cfg");
	}
	/* check file exists */
	fp = fopen(szCfgFile, "rb");
	if(fp)
	{
		fclose(fp);
		return 1;
	}
	return 0;
}

/* (slow) routines to find values from a cfg file */
int findKey(char*szCfgFile, char *szKey, char *szValue)
{
	int	 bRet;
	FILE *fp;

	bRet = 0;
	fp = fopen(szCfgFile,"rb");
	if(fp)
	{
		char szLine[N_MAXSTRLEN];

		while(!bRet && fgets(szLine,N_MAXSTRLEN,fp))
		{
			char *pPtr;
			if((pPtr = strchr(szLine,'='))!=NULL)
			{
				*pPtr = '\0';
				if(strcmp(szKey,szLine)==0)
				{
					strcpy(szValue,pPtr+1);
					/* strip off newline (if any) */
					if((pPtr = strpbrk(szValue,"\r\n"))!=NULL)
						*pPtr = '\0';
					bRet = 1;
				}
			}
		}
		fclose(fp);

	}
	return bRet;
}

int readStringFromCfg(char *szCfgFile, char *szKey, char *szValue)
{
    return(findKey(szCfgFile,szKey,szValue));
}

int readDoubleFromCfg(char *szCfgFile, char *szKey, double *pdValue)
{
	char szValue[N_MAXSTRLEN];

	if(findKey(szCfgFile,szKey,szValue))
	{
		*pdValue = atof(szValue);
		return 1;
	}
	return 0;
}

int readIntFromCfg(char *szCfgFile, char *szKey, int *pnValue)
{
	char szValue[N_MAXSTRLEN];

	if(findKey(szCfgFile,szKey,szValue))
	{
		*pnValue = atoi(szValue);
		return 1;
	}
	return 0;
}

