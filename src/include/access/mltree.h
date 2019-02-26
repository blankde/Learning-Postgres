/*
 * mltree.h
 *
 *  Created on: Sep 3, 2018
 *      Author: gpadmin
 */

#ifndef SRC_INCLUDE_ACCESS_MLTREE_H_
#define SRC_INCLUDE_ACCESS_MLTREE_H_


#include "access/genam.h"
#include "access/itup.h"
#include "access/sdir.h"
#include "access/xlog.h"
#include "access/xlogutils.h"

#include "access/nbtree.h"

//uint16 blkno = 0; //current page number
//uint16 blksize = 0;//total tuple of up to current page
//uint16* pageinfo;

/*
* mlpage like follows:
* +----------------+---------------------------------+
* | PageHeaderData | model1 | model2 | ...  | modelN |
* +-----------+----+---------------------------------+

* cause everyitem's size is same, we don't need postgres' Page Structure
* instead we use new mlpage structure
*/
typedef Page MLPage;

typedef PageHeaderData MLPageHeaderData;

typedef MLPageHeaderData *MLPageHeader;

/*
 * page-level and high-level locking modes (see README)
 */
#define ML_READ		BUFFER_LOCK_SHARE
#define ML_WRITE		BUFFER_LOCK_EXCLUSIVE
#define ML_NOLOCK		(-1)
#define ML_SHARE		ShareLock
#define ML_EXCLUSIVE	ExclusiveLock

#define MLEqualStrategyNumber			3

#define MAXMETABLOCK 512

/*
 *	We need to be able to tell the difference between read and write
 *	requests for pages, in order to do locking correctly.
 */

#define ML_READ			BUFFER_LOCK_SHARE
#define ML_WRITE		BUFFER_LOCK_EXCLUSIVE


#define MLMaxIndexTuplesPerPage \
		((int)(BLCKSZ - sizeof(MLPageHeaderData))/sizeof(MLModel))
//#define MLPageGetFreeSpace(page) \
	(BLCKSZ - sizeof(MLPageHeaderData) - ((MLPageHeader)page)->modelcount * sizeof(MLModel))
#define MLScanPosIsValid(scanpos) BufferIsValid((scanpos).buf)
#define MLPageGetMeta(page) \
	((MLMetaPage) PageGetContents(page))

typedef struct PageInfo{
	double blockmax;
	double blockmin;
	BlockNumber blkno;
} PageInfo;
typedef struct MLMetaPageData{
	uint16 length;
	bool isLastMeta;
	PageInfo pageInfo[MAXMETABLOCK];
}MLMetaPageData;

typedef MLMetaPageData *MLMetaPage;


/*
 * store data for ML model
 */
typedef struct MLModel{
	double t1;// sumOfXX = 0.0,
	double t2;// sumOfX = 0.0,
	double t3;// sumOfXY = 0.0
	double t4;// sumOfY = 0.0;
	double t5;// sumOfYY = 0.0;
	double k;
	double b;
	double maxerror;
	double minerror;
	double minvalue;
	double maxvalue;
	BlockNumber blkno;
	int length;
} MLModel;

typedef struct MLStateData{
	double xList[100000];
	double yList[100000];
	uint16 length;

	double rmse;
} MLStateData;


typedef uint16 SizeEntry;

typedef struct MLBuildState
{
	//bool		isUnique;
	//bool		haveDead;
	Relation	heapRel;
	BTSpool    *spool;//store (key,addr)

	/*
	 * spool2 is needed only when the index is an unique index. Dead tuples
	 * are put into spool2 instead of spool in order to avoid uniqueness
	 * check.
	 */
	//BTSpool    *spool2;
	double		indtuples;
	MLStateData* datastate;
	BlockNumber blkno;
	MLModel* mlmodel;

	uint16 modelsize;
	Datum xList[10000];
	bool* isValuesNull;
	ItemPointerData yList[100000];
	uint16 dataLength;

	BlockNumber pages_alloced;
	BlockNumber metablkno;

	Relation index;
	Relation	btreeIndex;

	Datum lastX; //we use this to record which values is biggest in precious block;
	Datum maxX; //we use this to record which values is biggest in current block;
} MLBuildState;

typedef struct MLPageState
{
	MLPage		mlps_page;		/* workspace for page building */
	BlockNumber mlps_blkno;		/* block # to write this page at */
	//IndexTuple	mlps_minkey;	/* copy of minimum key (first item) on page */
	//OffsetNumber mlps_lastoff;	/* last item offset loaded */
	//Size		mlps_full;		/* "full" if less than this much free space */
} MLPageState;

/*
 * Overall status record for index writing phase.
 */
typedef struct MLWriteState
{
	Relation index;
	MLMetaPage metap;
	Buffer metabuf;
	MLModel* mlmodel;
	uint16 modelsize;
	bool		mlws_use_wal;	/* dump pages to WAL? */
	BlockNumber mlws_pages_alloced;		/* # pages allocated */
	BlockNumber metablkno;
	//BlockNumber mlws_pages_written;		/* # pages written out, it's needed for non-sequencial writed. But now it's always sequcncial*/
	//Page		mlws_zeropage;	/* workspace for filling zeroes */
} MLWriteState;


typedef struct MLScanPosItem	/* what we remember about each match */
{
	ItemPointerData heapTid;	/* TID of referenced heap item */
	OffsetNumber indexOffset;	/* index item's location within page */
} MLScanPosItem;
typedef struct MLScanPosData
{
	Buffer		buf;			/* if valid, the buffer is pinned */

	BlockNumber nextPage;		/* page's right link when we scanned it */

	/*
	 * moreLeft and moreRight track whether we think there may be matching
	 * index entries to the left and right of the current page, respectively.
	 * We can clear the appropriate one of these flags when _bt_checkkeys()
	 * returns continuescan = false.
	 */
	bool		moreLeft;
	bool		moreRight;

	/*
	 * The items array is always ordered in index order (ie, increasing
	 * indexoffset).  When scanning backwards it is convenient to fill the
	 * array back-to-front, so we start at the last slot and fill downwards.
	 * Hence we need both a first-valid-entry and a last-valid-entry counter.
	 * itemIndex is a cursor showing which entry was last returned to caller.
	 */
	int			firstItem;		/* first valid index in items[] */
	int			lastItem;		/* last valid index in items[] */
	int			itemIndex;		/* current index in items[] */

	MLScanPosItem items[MLMaxIndexTuplesPerPage]; /* MUST BE LAST */
} MLScanPosData;

typedef struct MLScanOpaqueData
{
	/* these fields are set by _bt_preprocess_keys(): */
	bool		qual_ok;		/* false if qual can never be satisfied */
	int			numberOfKeys;	/* number of preprocessed scan keys */
	ScanKey		keyData;		/* array of preprocessed scan keys */

	/* info about killed items if any (killedItems is NULL if never used) */
	int		   *killedItems;	/* currPos.items indexes of killed items */
	int			numKilled;		/* number of currently stored items */

	/*
	 * If the marked position is on the same page as current position, we
	 * don't use markPos, but just keep the marked itemIndex in markItemIndex
	 * (all the rest of currPos is valid for the mark position). Hence, to
	 * determine if there is a mark, first look at markItemIndex, then at
	 * markPos.
	 */
	int			markItemIndex;	/* itemIndex, or -1 if not valid */

	/* keep these last in struct for efficiency */
	MLScanPosData currPos;		/* current position data */
	MLScanPosData markPos;		/* marked position, if any */
	bool isFisrtScaned; //
} MLScanOpaqueData;

typedef MLScanOpaqueData *MLScanOpaque;

/*
 * prototypes for functions in nbtree.c (external entry points for btree)
 */
extern Datum btbuild(PG_FUNCTION_ARGS);

//extern uint16 blkno; //current page number
//extern uint16 blksize;//total tuple of up to current page
//extern uint16* pageinfo;

#define MLPageGetModelCount(page) \
	((int)((PageHeader) (page)->pd_upper - sizeof(MLPageHeaderData))/sizeof(MLModel))

#define MLPageGetLastModel(page) \
	((MLModel*)((char *) (page)+((PageHeader) (page))->pd_lower - sizeof(MLModel)))

#define MLPageGetFirstModel(page) \
	((MLModel*)((char *) (page) + MAXALIGN(SizeOfPageHeaderData)))

#define MLPageGetModel(page,offset) \
	((MLModel*)((char *) (page) + MAXALIGN(SizeOfPageHeaderData)+offset*sizeof(MLModel)))

#endif /* SRC_INCLUDE_ACCESS_MLTREE_H_ */
