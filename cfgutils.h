#ifndef _CFGUTILS_H_
#define _CFGUTILS_H_

#define		N_MAXFNAMELEN			256				/* Maximum length of a filename */
#define		N_MAXSTRLEN				1024			/* Maximum length of a string */
#define		N_LINEBLOCKSIZE			256				/* How many lines to allocate at once */
#define		N_MAXREADINLEN			8192			/* Max length of input line */
/***********************/
/* platform specific   */
/***********************/
#ifdef _WIN32
#define 	C_DIR_DELIMITER '\\'
#else
#define 	C_DIR_DELIMITER '/'
#endif

int		getCfgFileName(char *szProgName, char *szCfgFile);
int 	readStringFromCfg(char *szCfgFile, char *szKey, char *szValue);
int 	readDoubleFromCfg(char *szCfgFile, char *szKey, double *pdValue);
int 	readIntFromCfg(char *szCfgFile, char *szKey, int *pnValue);

#endif /* _CFGUTILS_H_ */

