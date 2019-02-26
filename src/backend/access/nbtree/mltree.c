/*
 * mltree.c
 *
 *  Created on: Jul 4, 2018
 *      Author: gpadmin
 */

#include "postgres.h"

#include "access/genam.h"
#include "catalog/pg_am.h"
#include "access/mltree.h"
#include "access/relscan.h"
#include "catalog/index.h"
#include "catalog/storage.h"
#include "commands/vacuum.h"
#include "storage/bufmgr.h"
#include "storage/freespace.h"
#include "storage/indexfsm.h"
#include "storage/ipc.h"
#include "storage/lmgr.h"
#include "storage/smgr.h"
#include "tcop/tcopprot.h"
#include "utils/memutils.h"
#include "pgstat.h"
#include "storage/bufpage.h"
#include "access/heapam.h"
#include "catalog/indexing.h"
#include "catalog/pg_mlbtree.h"
#include "math.h"
#include "utils/fmgroids.h"
#include "utils/tqual.h"
//create index i1 on orders using mltree(o_orderkey);



	instr_time	before,
						after;
extern	double time1 = 0.0;
void _ls_insert(Datum value, uint16 offset,MLModel* model,Relation heapRel, Relation index){
	double this_error;

	this_error = model->k * value + model->b - offset;
	if(this_error > model->maxerror || this_error < model->minerror){
		double t1,t2,t3,t4,t5;
		Datum indexKey;
		Datum ip_posid;
		int2 length;
		double current_error;
		double max_error = 0;
		double min_error = 0;
		Datum xList[1000];
		Datum yList[1000];
		bool isnull;
		int2 keycol;



		t1 = t2 = t3 = t4 = t5 = 0.0;
		length = 0;
		Assert(index->rd_index->indnatts==1);
		keycol = index->rd_index->indkey.values[0];
	    heap_open(heapRel->rd_id,ML_READ);
		isnull = false;

		HeapTupleData	ctup;
		Buffer cbuf;
		Page		dp;
		int			lines;
		int			ntup;
		OffsetNumber lineoff;
		ItemId		lpp;


		ctup.t_tableOid = RelationGetRelid(heapRel);
		cbuf = ReadBufferExtended(heapRel, MAIN_FORKNUM, model->blkno,
											   RBM_NORMAL, NULL);
		dp = (Page) BufferGetPage(cbuf);
		lines = PageGetMaxOffsetNumber(dp);


		for (lineoff = FirstOffsetNumber, lpp = PageGetItemId(dp, lineoff);
				 lineoff <= lines;
				 lineoff++, lpp++){
			ctup.t_data = (HeapTupleHeader) PageGetItem((Page) dp, lpp);
			ctup.t_len = ItemIdGetLength(lpp);
			ItemPointerSet(&(ctup.t_self), model->blkno, lineoff);
			indexKey = heap_getattr(&ctup, keycol, RelationGetDescr(heapRel), &isnull);
			ip_posid = ctup.t_self.ip_posid;
			t1 += indexKey * indexKey;
			t2 += indexKey;
			t3 += indexKey * ip_posid;
			t4 += ip_posid;
			t5 += ip_posid * ip_posid;
			xList[length] = indexKey;
			yList[length] = ip_posid;
			length = length + 1;
		}
		ReleaseBuffer(cbuf);
		model->k = (t3*length - t2*t4) / (t1*length - t2*t2);
		model->b = (t4 - model->k* t2) / length;

		int i;
		for (i = 0; i<length;i++) {
			current_error =  yList[i] - model->k * xList[i] - model->b;
			max_error = (max_error > current_error) ? max_error:current_error;
			min_error = (min_error < current_error) ? min_error:current_error;
		}

		double error_mean,tempV,error_sigma;
		double k = model->k;
		double b = model->b;
		error_mean = (t4-k*t2)/length - b;
		tempV = error_mean + b;
		error_sigma = sqrt((t5+k*k*t1+2*k*tempV*t2-2*k*t3-2*tempV*t4)/length+tempV*tempV);
		model->maxerror = max_error+error_sigma;
		model->minerror = min_error-error_sigma;
		heap_close(heapRel,NoLock);
		//pfree(ctup);
	}

	model->maxvalue = (model->maxvalue > value) ? model->maxvalue:value;
	model->length++;
}


_init_model(MLModel* mlmodel){
	mlmodel->t1 = 0;
	mlmodel->t2 = 0;
	mlmodel->t3 = 0;
	mlmodel->t4 = 0;
	mlmodel->t5 = 0;
	mlmodel->k = 0;
	mlmodel->b = 0;
	mlmodel->minerror = 0;
	mlmodel->maxerror = 0;
}

void _ls_train(MLBuildState* state){
	MLModel* mlmodel =  &(state->mlmodel[state->modelsize++]);
	Datum* xList = state->xList;
	ItemPointerData* yList  = state->yList;
	uint16 length = state->dataLength;
	BlockNumber blkno = state->blkno;
	if(length == 0) return;
	double current_error;
	double max_error = 0;
	double min_error = 0;
	double deno;
	double error_sigma,error_mean, tempV;
	double t1,t2,t3,t4,t5;
	double k,b;
	int i;
	int count = length;
	double maxX = 0.0;
	double minX = 0.0;
	bool isunvalid[1000];

	t1 = t2 = t3 = t4 = t5 = 0.0;
	_init_model(mlmodel);

	if(length != 1){
		for(i = 0;i<length;i++){
			t1 += xList[i] * xList[i];
			t2 += xList[i];
			t3 += xList[i] * yList[i].ip_posid;
			t4 += yList[i].ip_posid;
			t5 += yList[i].ip_posid * yList[i].ip_posid;
		}
		if((deno = t1*length - t2*t2)==0) return;
		k = (t3*length - t2*t4) / deno;
		b = (t4 - k * t2) / length;
		error_mean = (t4-k*t2)/length - b;
		tempV = error_mean + b;
		error_sigma = sqrt((t5+k*k*t1+2*k*tempV*t2-2*k*t3-2*tempV*t4)/length+tempV*tempV);

		for (i = 0; i<count;i++) {
			current_error =  yList[i].ip_posid - k * xList[i] - b;
			if(current_error < error_mean - 3*error_sigma || current_error > error_mean + 3*error_sigma){
				IndexTuple	itup;

				t1 -= xList[i] * xList[i];
				t2 -= xList[i];
				t3 -= xList[i] * yList[i].ip_posid;
				t4 -= yList[i].ip_posid;
				t5 -= yList[i].ip_posid * yList[i].ip_posid;
				isunvalid[i] = true;
				length--;
				/* form an index tuple and point it at the heap tuple */
				itup = index_form_tuple(RelationGetDescr(state->btreeIndex), &xList[i], state->isValuesNull);
				itup->t_tid = yList[i];

				/*
				 * insert the index tuple into the appropriate spool file for subsequent
				 * processing
				 */
				_bt_spool(itup, state->spool);
			}
		}
		if((deno = t1*length - t2*t2)==0) return;
		mlmodel->k = (t3*length - t2*t4) / deno;
		mlmodel->b = (t4 - k * t2) / length;
		for (i = 0; i<count;i++) {
			if(isunvalid[i]) continue;
			maxX = (maxX>xList[i])?maxX:xList[i];
			minX = (minX<xList[i])?minX:xList[i];
			current_error =  yList[i].ip_posid - mlmodel->k * xList[i] - mlmodel->b;
			max_error = (max_error > current_error) ? max_error:current_error;
			min_error = (min_error < current_error) ? min_error:current_error;
		}
	}

	else{
		mlmodel->k = 1;
		mlmodel->b = yList[0].ip_posid - xList[0];
	}
	error_mean = (t4-k*t2)/length - b;
	tempV = error_mean + b;
	error_sigma = sqrt((t5+k*k*t1+2*k*tempV*t2-2*k*t3-2*tempV*t4)/length+tempV*tempV);
	mlmodel->t1 = t1;
	mlmodel->t2 = t2;
	mlmodel->t3 = t3;
	mlmodel->t4 = t4;
	mlmodel->t5 = t5;
	mlmodel->maxerror = max_error+error_sigma;
	mlmodel->minerror = min_error-error_sigma;
	mlmodel->maxvalue = state->maxX = maxX;
	mlmodel->minvalue =  minX;
	mlmodel->blkno = blkno;
	mlmodel->length = length;
}

void _ml_pushmodel(MLBuildState* state ){
	/*
	MLStateData* mldatastate = state->datastate;
	MLModel* mlmodel = state->mlmodel;

	double max_error = 0;
	double min_error = 0;
	double current_error;
	int i;
	//我们不需要计算最后一个数，因为最后一个数不再模型之中
	for (i = 0; i < mldatastate->length-1; i++) {
		current_error =  mldatastate->pre_factor[1][1] * mldatastate->xList[i] + mldatastate->pre_factor[0][1] - mldatastate->yList[i][1];
		max_error = (max_error > current_error) ? max_error:current_error;
		min_error = (min_error < current_error) ? min_error:current_error;
	}
	mlmodel[state->modelsize].k[0] = mldatastate->pre_factor[1][0];
	mlmodel[state->modelsize].k[1] = mldatastate->pre_factor[1][1];
	mlmodel[state->modelsize].b[0] = mldatastate->pre_factor[0][0];
	mlmodel[state->modelsize].b[1] = mldatastate->pre_factor[0][1];
	mlmodel[state->modelsize].maxerror = max_error;
	mlmodel[state->modelsize].minerror = min_error;
	mlmodel[state->modelsize].maxvalue = mldatastate->xList[mldatastate->length-2];
	state->modelsize = state->modelsize+1;
	*/
}

void _ml_pushlastmodel(MLBuildState* state ){
	/*
	MLStateData* mldatastate = state->datastate;
	MLModel* mlmodel = state->mlmodel;

	double max_error = 0;
	double min_error = 0;
	double current_error;
	int i;
	//我们不需要计算最后一个数，因为最后一个数不再模型之中
	for (i = 0; i < mldatastate->length-1; i++) {
		current_error =  mldatastate->pre_factor[1][1] * mldatastate->xList[i] + mldatastate->pre_factor[0][1] - mldatastate->yList[i][1];
		max_error = (max_error > current_error) ? max_error:current_error;
		min_error = (min_error < current_error) ? min_error:current_error;
	}
	mlmodel[state->modelsize].k[0] = mldatastate->pre_factor[1][0];
	mlmodel[state->modelsize].k[1] = mldatastate->pre_factor[1][1];
	mlmodel[state->modelsize].b[0] = mldatastate->pre_factor[0][0];
	mlmodel[state->modelsize].b[1] = mldatastate->pre_factor[0][1];
	mlmodel[state->modelsize].maxerror = max_error;
	mlmodel[state->modelsize].minerror = min_error;
	mlmodel[state->modelsize].maxvalue = mldatastate->xList[mldatastate->length-1];
	state->modelsize = state->modelsize+1;
	*/
}

static void mlbuildCallback(Relation index,
				HeapTuple htup,
				Datum *values,
				bool *isnull,
				bool tupleIsAlive,
				void *state);
/*
 * create and initialize a spool structure
 * maybe we need to sort it first,but now we don't do it,just like it's sorted.
 * The model should be palloed memory by the caller.
 */
MLBuildState *
_ml_init(Relation heap,Relation index, bool isunique, bool isdead)
{

	MLBuildState* state = (MLBuildState *) palloc0(sizeof(MLBuildState));

	instr_time	before,
						after;
	INSTR_TIME_SET_CURRENT(before);

	state->mlmodel= (MLModel *) palloc0((ALLOCSET_DEFAULT_INITSIZE)*sizeof(MLModel));
	state->modelsize= 0;


	INSTR_TIME_SET_CURRENT(after);
	INSTR_TIME_SUBTRACT(after, before);
	printf(_("Time: %.3f ms\n"), INSTR_TIME_GET_MILLISEC(after));


	state->dataLength = 0;

	state->blkno = -1;//Invalid
	state->pages_alloced = 0;
	state->indtuples = 0;
	state->index = index;

	return state;
}

void _ml_reset(MLBuildState *state){
	state->dataLength = 0;
}

/*
 *	_ml_pageinit() -- Initialize a new hash index page.
 */
void
_ml_pageinit(Page page, Size size)
{
	PageInit(page, size, 0);
}


MLPage _ml_newpage(){
	MLPage		page;
	MLPageHeader	p;

	page = (MLPage) palloc0(BLCKSZ);
	p = (MLPageHeader)page;

	/* Zero the page and set up standard page header info */
	_ml_pageinit(page, BLCKSZ);

	return page;
}


/*
 *	_ml_getnewbuf() -- Get a new page at the end of the index.
 *
 *		It is caller's responsibility to ensure that only one process can
 *		extend the index at a time.  In practice, this function is called
 *		only while holding write lock on the metapage, because adding a page
 *		is always associated with an update of metapage data.
 */
Buffer
_ml_getnewbuf(Relation rel, BlockNumber blkno, ForkNumber forkNum)
{
	BlockNumber nblocks = RelationGetNumberOfBlocksInFork(rel, forkNum);
	Buffer		buf;

	if (blkno == P_NEW)
		elog(ERROR, "ml AM does not use P_NEW");
	if (blkno > nblocks)
		elog(ERROR, "access to noncontiguous page in ml index \"%s\"",
			 RelationGetRelationName(rel));

	/* smgr insists we use P_NEW to extend the relation */
	if (blkno == nblocks)
	{
		buf = ReadBufferExtended(rel, forkNum, P_NEW, RBM_NORMAL, NULL);
		if (BufferGetBlockNumber(buf) != blkno)
			elog(ERROR, "unexpected ml relation size: %u, should be %u",
				 BufferGetBlockNumber(buf), blkno);
		LockBuffer(buf, ML_WRITE);
	}
	else
	{
		buf = ReadBufferExtended(rel, forkNum, blkno, RBM_ZERO_AND_LOCK,
								 NULL);
	}

	/* ref count and lock type are correct */

	/* initialize the page */
	_ml_pageinit(BufferGetPage(buf), BufferGetPageSize(buf));

	return buf;
}

/*
 *	_ml_wrtbuf() -- write a hash page to disk.
 *
 *		This routine releases the lock held on the buffer and our refcount
 *		for it.  It is an error to call _hash_wrtbuf() without a write lock
 *		and a pin on the buffer.
 *
 * NOTE: this routine should go away when/if hash indexes are WAL-ified.
 * The correct sequence of operations is to mark the buffer dirty, then
 * write the WAL record, then release the lock and pin; so marking dirty
 * can't be combined with releasing.
 */
void
_ml_wrtbuf(Relation rel, Buffer buf)
{
	MarkBufferDirty(buf);
	UnlockReleaseBuffer(buf);
}

/*
 *	_hash_metapinit() -- Initialize the metadata page of a hash index,
 *				the initial buckets, and the initial bitmap page.
 *
 * The initial number of buckets is dependent on num_tuples, an estimate
 * of the number of tuples to be loaded into the index initially.  The
 * chosen number of buckets is returned.
 *
 * We are fairly cavalier about locking here, since we know that no one else
 * could be accessing this index.  In particular the rule about not holding
 * multiple buffer locks is ignored.
 */

/* return the metabuf */
Buffer  _ml_metapinit(Relation index,BlockNumber blkno)
{
	Buffer metabuf;
	MLMetaPage metap;
	Page		pg;

	/*for now we don't deal ffactor*/

	metabuf = _ml_getnewbuf(index, blkno, MAIN_FORKNUM);
	pg = BufferGetPage(metabuf);


	metap = MLPageGetMeta(pg);
	metap->length = 0;
	metap->isLastMeta = true;

	MemSet(metap->pageInfo, 0, sizeof(metap->pageInfo));

	MarkBufferDirty(metabuf);

	return metabuf;
	/* all done */
	//_ml_wrtbuf(wstate->index, metabuf);
}


MLPageState* _ml_newpagestate(MLWriteState* wstate){
	MLPageState* state = (MLPageState *) palloc0(sizeof(MLPageState));
	state->mlps_page = _ml_newpage();
	state->mlps_blkno = wstate->mlws_pages_alloced++;

	return state;
}

void _ml_writepage(Relation index, MLPage page, BlockNumber blkno,ForkNumber forknum){
	//for now we don't consider transactions


	/*
	 * Now write the page.  Maybe it's always sequencial writing,so it's just need extend...I am not sure
	 */
	smgrextend(index->rd_smgr, forknum, blkno, (char *) page, true);
}

OffsetNumber MLPageAddItem(Page page,Item item,Size size){
	PageHeader	phdr = (PageHeader) page;
	memcpy((char *) page + phdr->pd_lower, item, size);
	phdr->pd_lower += size;
}
void ml_addtopage(MLWriteState* wstate, MLPageState* state,MLModel* mlmodel){

	MLPage npage = state->mlps_page;
	BlockNumber nblkno = state->mlps_blkno;
	Size modelsize = sizeof(MLModel);
	Size pgspc = PageGetFreeSpace(npage);

	 /* Check whether the model can fit on a mlree page at all.*/

	if(pgspc < modelsize){

		 /* Finish off the page and write it out.*/

		MLPage		opage = npage;
		BlockNumber oblkno = nblkno;
		npage = _ml_newpage();

		//page is up to MAXMETABLOCK,we alloc a new metapage to record meta info.
		if((nblkno+1) % MAXMETABLOCK==0){
			//maybe here we don't need to mark the buffer dirty because it's dirty
			wstate->metap->isLastMeta = false;
			UnlockReleaseBuffer(wstate->metabuf);
			wstate->metabuf = _ml_metapinit(wstate->index,wstate->mlws_pages_alloced);
			wstate->metap = MLPageGetMeta(BufferGetPage(wstate->metabuf));
			wstate->metablkno = wstate->mlws_pages_alloced++;
		}
		//move to next page
		nblkno = wstate->mlws_pages_alloced++;
		/* Write out the old page.  We never need to touch it again, so we can
		 * free the opage workspace too.*/
		wstate->metap->pageInfo[wstate->metap->length].blockmax = MLPageGetLastModel(state->mlps_page)->maxvalue;
		wstate->metap->pageInfo[wstate->metap->length].blockmin = MLPageGetLastModel(state->mlps_page)->minvalue;
		wstate->metap->pageInfo[wstate->metap->length++].blkno = MLPageGetLastModel(state->mlps_page)->blkno;
		MarkBufferDirty(wstate->metabuf);
		_ml_writepage(wstate->index, opage, oblkno,MAIN_FORKNUM);

		state->mlps_page = npage;
		state->mlps_blkno = nblkno;
	}

	 /* Add the new item into the current page.*/


	MLPageAddItem(npage,(Item)mlmodel,modelsize);
}

void ml_write(MLBuildState* buildstate){
	//Buffer		metabuf;
	//MLMetaPage		metapage;
	Buffer		buf;
	Page		page;
	MLPageState* state;
	int i;

	uint16 modelsize = buildstate->modelsize;

	MLWriteState wstate;
	wstate.mlmodel = buildstate->mlmodel;
	wstate.modelsize = buildstate->modelsize;
	wstate.mlws_pages_alloced = buildstate->pages_alloced;
	wstate.metablkno = buildstate->metablkno;
	//wstate.mlws_pages_written = 0;
	//wstate.mlws_zeropage = NULL;	/* until needed */
	wstate.index = buildstate->index;

	//alloc the first meta page
	if(wstate.mlws_pages_alloced == 0){
		wstate.metabuf = _ml_metapinit(wstate.index,0);
		wstate.metap = MLPageGetMeta(BufferGetPage(wstate.metabuf));
		wstate.metablkno = wstate.mlws_pages_alloced++;

		state= _ml_newpagestate(&wstate);
	}
	else{
		wstate.metabuf = _ml_getbuf(buildstate->index, buildstate->metablkno, ML_READ);
		wstate.metap = MLPageGetMeta(BufferGetPage(wstate.metabuf));

		buf = _ml_getbuf(buildstate->index, wstate.metap->length, ML_READ);
		page = BufferGetPage(buf);

		if(PageGetFreeSpace(page) < sizeof(MLModel)){
			state= _ml_newpagestate(&wstate);
		}
		else{
			state = (MLPageState *) palloc0(sizeof(MLPageState));
			state->mlps_page = page;
			state->mlps_blkno = wstate.mlws_pages_alloced-1;
			wstate.metap->length--;
		}
	}
	//add model to page
	for(i = 0;i<modelsize;i++){
		ml_addtopage(&wstate,state,&wstate.mlmodel[i]);
	}

	//write last page
	MLModel* m = MLPageGetLastModel(state->mlps_page);
	wstate.metap->pageInfo[wstate.metap->length].blockmax = MLPageGetLastModel(state->mlps_page)->maxvalue;
	wstate.metap->pageInfo[wstate.metap->length].blockmin = MLPageGetLastModel(state->mlps_page)->minvalue;
	wstate.metap->pageInfo[wstate.metap->length++].blkno = MLPageGetLastModel(state->mlps_page)->blkno;
	wstate.metap->isLastMeta = true;
	_ml_wrtbuf(wstate.index, wstate.metabuf);
	_ml_writepage(wstate.index, state->mlps_page, state->mlps_blkno,MAIN_FORKNUM);
	if(BufferIsValid(buf))
		UnlockReleaseBuffer(buf);
	FlushRelationBuffers(wstate.index);

	buildstate->pages_alloced = wstate.mlws_pages_alloced;
	buildstate->metablkno = wstate.metablkno;

	pfree(state->mlps_page);
	pfree(state);

	//now write the ipc file
	//FileWrite(wstate->rd_smgr, INDEX_PAGE_COUNT, blkno,(char *) pageinfo, true);
}

Relation
heap_create1(const char *relname,
			Oid relnamespace,
			Oid reltablespace,
			Oid relid,
			TupleDesc tupDesc,
			char relkind,
			char relpersistence,
			bool shared_relation,
			bool mapped_relation,
			bool allow_system_table_mods)
{
	bool		create_storage;
	Relation	rel;

	/* The caller must have provided an OID for the relation. */
	Assert(OidIsValid(relid));

	/*
	 * Decide if we need storage or not, and handle a couple other special
	 * cases for particular relkinds.
	 */
	switch (relkind)
	{
		case RELKIND_VIEW:
		case RELKIND_COMPOSITE_TYPE:
		case RELKIND_FOREIGN_TABLE:
			create_storage = false;

			/*
			 * Force reltablespace to zero if the relation has no physical
			 * storage.  This is mainly just for cleanliness' sake.
			 */
			reltablespace = InvalidOid;
			break;
		case RELKIND_SEQUENCE:
			create_storage = true;

			/*
			 * Force reltablespace to zero for sequences, since we don't
			 * support moving them around into different tablespaces.
			 */
			reltablespace = InvalidOid;
			break;
		default:
			create_storage = true;
			break;
	}

	/*
	 * Never allow a pg_class entry to explicitly specify the database's
	 * default tablespace in reltablespace; force it to zero instead. This
	 * ensures that if the database is cloned with a different default
	 * tablespace, the pg_class entry will still match where CREATE DATABASE
	 * will put the physically copied relation.
	 *
	 * Yes, this is a bit of a hack.
	 */
	if (reltablespace == 0)
		reltablespace = InvalidOid;

	/*
	 * build the relcache entry.
	 */

	rel = RelationBuildLocalRelation(relname,
									 relnamespace,
									 tupDesc,
									 relid,
									 reltablespace,
									 shared_relation,
									 mapped_relation,
									 relpersistence);

	/*
	 * Have the storage manager create the relation's disk file, if needed.
	 *
	 * We only create the main fork here, other forks will be created on
	 * demand.
	 */

	if (create_storage)
	{
		RelationOpenSmgr(rel);
		RelationCreateStorage(rel->rd_node, relpersistence);
	}


	return rel;
}

Datum
mlbuild(PG_FUNCTION_ARGS)
{
	Relation	heap = (Relation) PG_GETARG_POINTER(0);
	Relation	index = (Relation) PG_GETARG_POINTER(1);
	IndexInfo  *indexInfo = (IndexInfo *) PG_GETARG_POINTER(2);
	IndexBuildResult *result;
	double		reltuples;
	MLBuildState* buildstate;

	Relation pg_class;
	Relation pg_mlbtree;
	Relation pg_attribute;
	CatalogIndexState indstate;
	Oid indexRelationId;
	Oid tableSpaceId;
	char relpersistence;
	TupleDesc	indexTupDesc;
	TupleDesc	btreeTupDesc;

	INSTR_TIME_SET_CURRENT(before);
	/*
	 * Create the file first if it doesn't exist.  If smgr_vm_nblocks is
	 * positive then it must exist, no need for an smgrexists call.
	 */
	//if (!smgrexists(index->rd_smgr, INDEX_PAGE_COUNT))
		//smgrcreate(index->rd_smgr, INDEX_PAGE_COUNT, false);
	/*
	 * We expect to be called exactly once for any index relation. If that's
	 * not the case, big trouble's what we have.
	 */
	if (RelationGetNumberOfBlocks(index) != 0)
		elog(ERROR, "index \"%s\" already contains data",
			 RelationGetRelationName(index));

/*	//create the btree
	pg_class = heap_open(RelationRelationId, RowExclusiveLock);
	pg_mlbtree = heap_open(MLBTREERelationId, RowExclusiveLock);
	relpersistence = index->rd_rel->relpersistence;
	char* indexRelationName = strcat(index->rd_rel->relname.data,"_btree");
	tableSpaceId = index->rd_rel->reltablespace;
	indexTupDesc = CreateTupleDescCopy(index->rd_att);

	indexRelationId =
			GetNewRelFileNode(tableSpaceId, pg_class, relpersistence);

	btreeIndex = heap_create1(indexRelationName,
							index->rd_rel->relnamespace,
							tableSpaceId,
							indexRelationId,
							indexTupDesc,
							RELKIND_INDEX,
							relpersistence,
							false,
							false,
							false);
	Assert(indexRelationId == RelationGetRelid(btreeIndex));
	LockRelation(btreeIndex, AccessExclusiveLock);


	 * Fill in fields of the index's pg_class entry that are not set correctly
	 * by heap_create.
	 *
	 * XXX should have a cleaner way to create cataloged indexes


	btreeIndex->rd_rel->relowner = index->rd_rel->relowner;
	btreeIndex->rd_rel->relam = 403; //btree
	btreeIndex->rd_rel->relkind = RELKIND_INDEX;
	btreeIndex->rd_rel->relhasoids = false;


	 * store index's pg_class entry


	InsertPgClassTuple(pg_class, btreeIndex,
					   RelationGetRelid(btreeIndex),
					   (Datum) 0,
					   NULL);
	InsertPgMlBtreeTuple(pg_mlbtree,index->rd_id,btreeIndex->rd_id);

	//insert into pg_attribute
	pg_attribute = heap_open(AttributeRelationId, RowExclusiveLock);

	indstate = CatalogOpenIndexes(pg_attribute);


	 * insert data from new index's tupdesc into pg_attribute

	btreeTupDesc = RelationGetDescr(btreeIndex);

	for (int i = 0; i < index->rd_rel->relnatts; i++)
	{
		InsertPgAttributeTuple(pg_attribute, btreeTupDesc->attrs[i], indstate);
	}

	CatalogCloseIndexes(indstate);

	heap_close(pg_attribute, RowExclusiveLock);

	 done with pg_class

	heap_close(pg_class, RowExclusiveLock);*/




	buildstate = _ml_init(heap,index, indexInfo->ii_Unique, false);

	buildstate->btreeIndex = index_open(index->btreeIndexOid, RowExclusiveLock);

	buildstate->spool = _bt_spoolinit(buildstate->btreeIndex, true, false);


	/* do the heap scan */
		reltuples = IndexBuildHeapScan(heap, index, indexInfo, true,true,
									   mlbuildCallback, buildstate);


	//write the last model
	_ls_train(buildstate);
	ml_write(buildstate);

	result = (IndexBuildResult *) palloc(sizeof(IndexBuildResult));


	result->heap_tuples = reltuples;
	result->index_tuples = buildstate->modelsize;

	_bt_leafbuild(buildstate->spool, NULL);
	_bt_spooldestroy(buildstate->spool);
	index_close(buildstate->btreeIndex, NoLock);
	//pfree(buildstate->datastate);
	pfree(buildstate->mlmodel);
	pfree(buildstate);
	INSTR_TIME_SET_CURRENT(after);
	INSTR_TIME_SUBTRACT(after, before);
	time1 += INSTR_TIME_GET_MILLISEC(after);
	PG_RETURN_POINTER(result);
}

/*
 * Per-tuple callback from IndexBuildHeapScan
 *
 * only for 1d
 */
static void
mlbuildCallback(Relation index,
				HeapTuple htup,
				Datum *values,
				bool *isnull,
				bool tupleIsAlive,
				void *state)
{
	MLBuildState *buildstate = (MLBuildState *) state;
	IndexTuple	itup;
	uint16 bi_hi,bi_lo,pos;
	uint32 blkno;


	bi_hi = htup->t_self.ip_blkid.bi_hi;
	bi_lo = htup->t_self.ip_blkid.bi_lo;
	pos = htup->t_self.ip_posid;

	blkno = (bi_hi << 16) + bi_lo;
	if(buildstate->blkno == -1)
		buildstate->blkno = blkno;
	buildstate->isValuesNull = isnull;

	buildstate->indtuples += 1;

	if(buildstate->dataLength<10000 && blkno == buildstate->blkno){
		if(values[0]<buildstate->lastX){
			IndexTuple	itup;
			/* form an index tuple and point it at the heap tuple */
			itup = index_form_tuple(RelationGetDescr(buildstate->btreeIndex), values, isnull);
			itup->t_tid = htup->t_self;

			/*
			 * insert the index tuple into the appropriate spool file for subsequent
			 * processing
			 */
			_bt_spool(itup, buildstate->spool);
			return;
		}
		buildstate->xList[buildstate->dataLength] = values[0];
		buildstate->yList[buildstate->dataLength++] = htup->t_self;
		return;
	}

	//if blkno change, compute the model and save it
	_ls_train(buildstate);
	buildstate->dataLength = 0;
	buildstate->xList[buildstate->dataLength] = values[0];
	buildstate->yList[buildstate->dataLength++] = htup->t_self;
	buildstate->blkno = blkno;
	buildstate->lastX = buildstate->maxX;
	if(buildstate->modelsize >= ALLOCSET_DEFAULT_INITSIZE){
		ml_write(buildstate);
		buildstate->modelsize = 0;
	}
}


/*
 *	mlbeginscan() -- start a scan on a btree index
 */
Datum
mlbeginscan(PG_FUNCTION_ARGS)
{
	Relation	rel = (Relation) PG_GETARG_POINTER(0);
	int			nkeys = PG_GETARG_INT32(1);
	int			norderbys = PG_GETARG_INT32(2);
	IndexScanDesc scan;
	MLScanOpaque so;

	/* no order by operators allowed */
	Assert(norderbys == 0);

	/* get the scan */
	scan = RelationGetIndexScan(rel, nkeys, norderbys);

	/* allocate private workspace */
	so = (MLScanOpaque) palloc(sizeof(MLScanOpaqueData));
	so->currPos.buf = so->markPos.buf = InvalidBuffer;
	so->isFisrtScaned = true;
	if (scan->numberOfKeys > 0)
		so->keyData = (ScanKey) palloc(scan->numberOfKeys * sizeof(ScanKeyData));
	else
		so->keyData = NULL;
	so->killedItems = NULL;		/* until needed */
	so->numKilled = 0;
	scan->opaque = so;

	PG_RETURN_POINTER(scan);
}

/*
 *	mlrescan() -- rescan an index relation
 */
Datum
mlrescan(PG_FUNCTION_ARGS)
{
	IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
	ScanKey		scankey = (ScanKey) PG_GETARG_POINTER(1);

	/* remaining arguments are ignored */
	MLScanOpaque so = (MLScanOpaque) scan->opaque;

	/* we aren't holding any read locks, but gotta drop the pins */
	if (MLScanPosIsValid(so->currPos))
	{
		/* Before leaving current page, maybe we should deal with any killed items but now it doen't*/
		ReleaseBuffer(so->currPos.buf);
		so->currPos.buf = InvalidBuffer;
	}

	if (MLScanPosIsValid(so->markPos))
	{
		ReleaseBuffer(so->markPos.buf);
		so->markPos.buf = InvalidBuffer;
	}
	so->markItemIndex = -1;

	/*
	 * Reset the scan keys. Note that keys ordering stuff moved to _ml_first.
	 * - vadim 05/05/97
	 */
	if (scankey && scan->numberOfKeys > 0)
		memmove(scan->keyData,
				scankey,
				scan->numberOfKeys * sizeof(ScanKeyData));
	so->numberOfKeys = 0;		/* until _ml_preprocess_keys sets it */

	PG_RETURN_VOID();
}

/*
 *	mlendscan() -- close down a scan
 */
Datum
mlendscan(PG_FUNCTION_ARGS)
{
	IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
	MLScanOpaque so = (MLScanOpaque) scan->opaque;

	/* we aren't holding any read locks, but gotta drop the pins */
	if (MLScanPosIsValid(so->currPos))
	{
		/* Before leaving current page, maybe we should deal with any killed items but now it doen't */
		ReleaseBuffer(so->currPos.buf);
		so->currPos.buf = InvalidBuffer;
	}

	if (MLScanPosIsValid(so->markPos))
	{
		ReleaseBuffer(so->markPos.buf);
		so->markPos.buf = InvalidBuffer;
	}
	so->markItemIndex = -1;

	if (so->killedItems != NULL)
		pfree(so->killedItems);
	if (so->keyData != NULL)
		pfree(so->keyData);
	pfree(so);

	PG_RETURN_VOID();
}


/*
 *	_ml_getbuf() -- Get a buffer by block number for read or write.
 *
 *		'access' must be ML_READ, ML_WRITE, or ML_NOLOCK.
 *		'flags' is a bitwise OR of the allowed page types.
 *
 *
 *		When this routine returns, the appropriate lock is set on the
 *		requested buffer and its reference count has been incremented
 *		(ie, the buffer is "locked and pinned").
 *
 *		P_NEW is disallowed because this routine can only be used
 *		to access pages that are known to be before the filesystem EOF.
 *		Extending the index should be done with _ml_getnewbuf.
 */
Buffer
_ml_getbuf(Relation rel, BlockNumber blkno, int access)
{
	Buffer		buf;

	if (blkno == P_NEW)
		elog(ERROR, "ML AM does not use P_NEW");

	buf = ReadBuffer(rel, blkno);
	if (access != ML_NOLOCK)
		LockBuffer(buf, access);

	return buf;
}

int _ml_binarysearch(PageInfo arr[], double num, int sz)//二分查找函数
{
	if(num<arr[0].blockmax)
		return 0;
    int mid = 0;
    int left = 0;
    int right = sz - 1;
    while (left <= right)
    {
        mid = (left + right)>>1;
        if (num <arr[mid].blockmax)
            right = mid-1;
        else if (num > arr[mid].blockmax)
            left = mid+1;
        else
            return mid;
    }
    return (arr[mid].blockmax>=num) ? mid:mid+1;
}

int MLPageGetModelOffset(void* address,MLPage page){
	Assert(address>(void*)MLPageGetFirstModel(page));
	Assert(address<(void*)MLPageGetLastModel(page));

	return (address-(void*)MLPageGetFirstModel(page))/sizeof(MLModel);

}

/*
 *	btgettuple() -- Get the next tuple in the scan.
 */
Datum
btmlgettuple(IndexScanDesc scan,ScanDirection dir,BTScanOpaque so)
{
	bool		res;


	/* btree indexes are never lossy */
	scan->xs_recheck = false;

	/*
	 * If we've already initialized this scan, we can just advance it in the
	 * appropriate direction.  If we haven't done so yet, we call a routine to
	 * get the first item in the scan.
	 */
	if (BTScanPosIsValid(so->currPos))
	{
		/*
		 * Check to see if we should kill the previously-fetched tuple.
		 */
		if (scan->kill_prior_tuple)
		{
			/*
			 * Yes, remember it for later.  (We'll deal with all such tuples
			 * at once right before leaving the index page.)  The test for
			 * numKilled overrun is not just paranoia: if the caller reverses
			 * direction in the indexscan then the same item might get entered
			 * multiple times.  It's not worth trying to optimize that, so we
			 * don't detect it, but instead just forget any excess entries.
			 */
			if (so->killedItems == NULL)
				so->killedItems = (int *)
					palloc(MaxIndexTuplesPerPage * sizeof(int));
			if (so->numKilled < MaxIndexTuplesPerPage)
				so->killedItems[so->numKilled++] = so->currPos.itemIndex;
		}

		/*
		 * Now continue the scan.
		 */
		res = _bt_next(scan, dir);
	}
	else
		res = _bt_first(scan, dir);


	PG_RETURN_BOOL(res);
}


/*
 *	_ml_first() -- Find the first item in a scan.
 *
 *		Find the first item in the index that
 *		satisfies the qualification associated with the scan descriptor. On
 *		success, the page containing the current index tuple is read locked
 *		and pinned, and the scan's opaque data entry is updated to
 *		include the buffer.
 */
bool
_ml_first(IndexScanDesc scan, ScanDirection dir,BlockNumber blkno)
{
	Relation	rel = scan->indexRelation;
	MLScanOpaque so = (MLScanOpaque) scan->opaque;

	Buffer		buf;
	MLPage		page;
	ItemPointer current;
	ScanKey		cur;
	BlockNumber vblkno; //the number of block which store the predicated data
	int2 keycol;
	MLModel* targetModel;
	bool found = false;

	pgstat_count_index_scan(rel);

	current = &(so->currPos);

	/* There may be more than one index qual, but we only deal the first */
	cur = &scan->keyData[0];
	/* We support only single-column hash indexes */
	Assert(cur->sk_attno == 1);
	/* And there's only one operator strategy, too */
	Assert(cur->sk_strategy == MLEqualStrategyNumber);
	/*
	 * If the constant in the index qual is NULL, assume it cannot match any
	 * items in the index.
	 */
	if (cur->sk_flags & SK_ISNULL)
		return false;

	//find the targeted page which contains the model that we need.


	/*
	 * Acquire shared split lock
	 */
	if ((!RELATION_IS_LOCAL(rel)))
		LockPage(rel, blkno, ML_SHARE);
	/* Fetch the page*/
	buf = _ml_getbuf(rel, blkno, ML_READ);
	page = BufferGetPage(buf);
	targetModel = MLPageGetFirstModel(page);


	//do binary search on targeted page
    MLModel* left = MLPageGetFirstModel(page);
    MLModel* right = MLPageGetLastModel(page);
    MLModel* mid;
    if(cur->sk_argument>left->maxvalue){
    	while (left <= right)
    	{
    		mid = left + ((right - left)>>1);
    		if (cur->sk_argument <mid->maxvalue)
    			right = mid-1;
    		else if (cur->sk_argument > mid->maxvalue)
    			left = mid+1;
    		else{
    			break;
    		}
    	}

    	targetModel = (mid->maxvalue>=cur->sk_argument) ? mid:mid+1;
    }
    //HeapScanDesc heapscan;
    //HeapTuple	heapTuple;
   // TupleTableSlot *slot;
    //bool		isNull;
    //Datum		iDatum;
	//Snapshot snapshot;
    //ok,we find the model
    vblkno = targetModel->blkno;

	//get the predicated block
	Assert(rel->rd_index->indnatts==1);
	keycol = rel->rd_index->indkey.values[0];
    //heap_open(scan->heapRelation->rd_id,ML_READ);
/*
	snapshot = RegisterSnapshot(GetTransactionSnapshot());
	heapscan = heap_beginscan(scan->heapRelation,snapshot,1,cur);
	heapscan->rs_startblock = vblkno;
	heapscan->rs_cblock = vblkno;
	heapscan->rs_nblocks = vblkno+1;
	slot = MakeSingleTupleTableSlot(RelationGetDescr(scan->heapRelation)); //!!!
*/

	//now we do sequencial scan
	HeapTupleData	ctup;
	Buffer cbuf;
	Page		dp;
	int			lines;
	int			ntup;
	bool isnull;
	OffsetNumber lineoff;
	ItemId		lpp;
	OffsetNumber upper,lower;
	double predited;

	isnull = false;
	ctup.t_tableOid = RelationGetRelid(scan->heapRelation);
	cbuf = ReadBufferExtended(scan->heapRelation, MAIN_FORKNUM, vblkno,
										   RBM_NORMAL, NULL);
	dp = (Page) BufferGetPage(cbuf);
	lines = PageGetMaxOffsetNumber(dp);

	predited = targetModel->k * cur->sk_argument + targetModel->b;
	upper = floor(predited + targetModel->maxerror);
	lower = ceil(predited + targetModel->minerror);
	upper = upper>lines?lines:upper;
	lower = lower<0?0:lower;

	for (lineoff = lower, lpp = PageGetItemId(dp, lineoff);
			 lineoff <= upper;
			 lineoff++, lpp++){
		ctup.t_data = (HeapTupleHeader) PageGetItem((Page) dp, lpp);
		ctup.t_len = ItemIdGetLength(lpp);
		ItemPointerSet(&(ctup.t_self), vblkno, lineoff);
		if(cur->sk_argument == heap_getattr(&ctup, keycol, RelationGetDescr(scan->heapRelation), &isnull)){
			scan->xs_ctup.t_self = ctup.t_self;
			found = true;
			break;
		}
	}
	//if((heapTuple = heap_getpagenext(heapscan, ForwardScanDirection)) != NULL){
	//	scan->xs_ctup.t_self = heapTuple->t_self;
	//	found = true;
	//}
	if(!found){
		HeapTuple	heapTuple;
		Relation pg_mlbtree;
		HeapScanDesc pg_mlbtreescan;
		TupleDesc pg_mlbtreeDesc;
		Oid btreeIndexOid;
		Relation btreeIndex;
		IndexScanDesc scan2;
		HeapTuple pg_mlbtreeTuple;
		ScanKeyData key;
		Snapshot snapshot;

		snapshot = RegisterSnapshot(GetTransactionSnapshot());
		pg_mlbtree = heap_open(MLBTREERelationId, RowExclusiveLock);
		ScanKeyInit(&key,
					Anum_pg_mlbtree_ml,
					BTEqualStrategyNumber, F_OIDEQ,
					ObjectIdGetDatum(scan->indexRelation->rd_id));

		/* see notes above about using SnapshotDirty */
		pg_mlbtreescan = systable_beginscan(pg_mlbtree, RelationGetOidIndex(pg_mlbtree), true,
								snapshot, 1, &key);
		pg_mlbtreeDesc = RelationGetDescr(pg_mlbtree);
		if((pg_mlbtreeTuple = systable_getnext(pg_mlbtreescan)) != NULL){
			btreeIndexOid = DatumGetObjectId(nocachegetattr(pg_mlbtreeTuple,2,pg_mlbtreeDesc));
		}
		//heap_endscan(pg_mlbtreescan);
		systable_endscan(pg_mlbtreescan);
		heap_close(pg_mlbtree, NoLock);
		btreeIndex = index_open(btreeIndexOid, RowExclusiveLock);

		scan2 = index_beginscan(scan->heapRelation,btreeIndex,snapshot,scan->numberOfKeys,scan->numberOfOrderBys);
		index_rescan(scan2,scan->keyData,scan->numberOfKeys,NULL,scan->numberOfOrderBys);
		heapTuple = index_getnext(scan2,dir);
		if(heapTuple!=NULL){
			scan->xs_ctup.t_self = heapTuple->t_self;
			found = true;
		}
		index_endscan(scan2);
		index_close(btreeIndex, RowExclusiveLock);
		UnregisterSnapshot(snapshot);
	}

	/* done with pg_class */
	////UnregisterSnapshot(snapshot);
	//ExecDropSingleTupleTableSlot(slot);
	//heap_endscan(heapscan);
	UnlockReleaseBuffer(buf);
	ReleaseBuffer(cbuf);
	UnlockPage(rel, blkno, ML_SHARE);
	return found;

}

/*
 *	mlgettuple() -- Get the next tuple in the scan.
 */
Datum
mlgettuple(PG_FUNCTION_ARGS)
{
	IndexScanDesc scan = (IndexScanDesc) PG_GETARG_POINTER(0);
	ScanDirection dir = (ScanDirection) PG_GETARG_INT32(1);
	MLScanOpaque so = (MLScanOpaque) scan->opaque;
	bool		res;
	Buffer 	metabuf;
	MLMetaPage metap;
	int i = 0;
	INSTR_TIME_SET_CURRENT(before);

	/* mlree indexes are never lossy */
	scan->xs_recheck = false;


	/*
	 * If we've already initialized this scan, we can just advance it in the
	 * appropriate direction.  If we haven't done so yet, we call a routine to
	 * get the first item in the scan.
	 */
	if (so->isFisrtScaned)
	{
		BlockNumber blkno;
		metabuf = _ml_getbuf(scan->indexRelation, 0, ML_READ);
		metap = MLPageGetMeta(BufferGetPage(metabuf));

		while((scan->keyData[0].sk_argument)>(metap->pageInfo[metap->length-1].blockmax)){
			if(metap->isLastMeta)
				break;
			i++;
			UnlockReleaseBuffer(metabuf);
			metabuf = _ml_getbuf(scan->indexRelation, i*MAXMETABLOCK, ML_READ);
			metap = MLPageGetMeta(BufferGetPage(metabuf));
		}
		blkno = (BlockNumber)_ml_binarysearch(metap->pageInfo,scan->keyData[0].sk_argument,metap->length);
		blkno = i*MAXMETABLOCK+blkno+1;
		so->isFisrtScaned = false;

		INSTR_TIME_SET_CURRENT(after);
		INSTR_TIME_SUBTRACT(after, before);

		time1 += INSTR_TIME_GET_MILLISEC(after);
		res = _ml_first(scan, dir,blkno);

		UnlockReleaseBuffer(metabuf);
	}
	else{
		res = false;
	}

	PG_RETURN_BOOL(res);
}

/* using binary search,return where the proper model of the block locates*/
BlockNumber _ml_findTargetPage(PageInfo arr[],BlockNumber blkno,BlockNumber totalCount) {
	if(blkno<arr[0].blkno)
		return 0;
    int mid = 0;
    int left = 0;
    int right = totalCount - 1;
    while (left <= right)
    {
        mid = (left + right)>>1;
        if (blkno <arr[mid].blkno)
            right = mid-1;
        else if (blkno > arr[mid].blkno)
            left = mid+1;
        else
            return mid;
    }
    return (arr[mid].blkno>=blkno) ? mid:mid+1;
}

Datum
mlinsert(PG_FUNCTION_ARGS)
{
	INSTR_TIME_SET_CURRENT(before);
	//INSTR_TIME_SET_CURRENT(before2);
	Relation	rel = (Relation) PG_GETARG_POINTER(0);
	Datum	   *values = (Datum *) PG_GETARG_POINTER(1);
	bool	   *isnull = (bool *) PG_GETARG_POINTER(2);
	ItemPointer ht_ctid = (ItemPointer) PG_GETARG_POINTER(3);
	Relation	heapRel = (Relation) PG_GETARG_POINTER(4);
	IndexUniqueCheck checkUnique = (IndexUniqueCheck) PG_GETARG_INT32(5);
	bool		result;
	Buffer 	metabuf;
	MLMetaPage metap;
	Buffer buf;
	Page page;
	MLModel* lastModel;
	BlockNumber metablkno,targetBlkno;
	MLModel* targetModel;
	int i = 0;

	uint16 bi_hi,bi_lo,pos;
	uint32 blkno;

	//INSTR_TIME_SET_CURRENT(before);
	bi_hi = ht_ctid->ip_blkid.bi_hi;
	bi_lo = ht_ctid->ip_blkid.bi_lo;
	blkno = (bi_hi << 16) + bi_lo;
	pos = ht_ctid->ip_posid;

	metabuf = _ml_getbuf(rel, 0, ML_READ);
	metap = MLPageGetMeta(BufferGetPage(metabuf));
	while(blkno>(metap->pageInfo[metap->length-1].blkno)){
		if(metap->isLastMeta)
			break;
		i++;
		UnlockReleaseBuffer(metabuf);
		metabuf = _ml_getbuf(rel, i*MAXMETABLOCK, ML_READ);
		metap = MLPageGetMeta(BufferGetPage(metabuf));
	}
	//INSTR_TIME_SET_CURRENT(after);
	//INSTR_TIME_SUBTRACT(after, before);
	//time4 += INSTR_TIME_GET_MILLISEC(after);
	if(blkno>metap->pageInfo[metap->length-1].blkno){
		//INSTR_TIME_SET_CURRENT(before);
		MLModel* mlmodel;
		MLPage		npage;
		BlockNumber oblkno, nblkno;

		mlmodel = palloc0(sizeof(MLModel));
		mlmodel->t1 = values[0]*values[0];
		mlmodel->t2 = values[0];
		mlmodel->t3 = values[0]*pos;
		mlmodel->t4 = pos;
		mlmodel->k = 1;
		mlmodel->b = pos - DatumGetFloat8(values[0]);
		mlmodel->maxvalue = values[0];
		mlmodel->blkno = blkno;
		mlmodel->length = 1;


		targetBlkno = i*MAXMETABLOCK+metap->length;
		/*
		 * Acquire shared split lock
		 */
		if ((!RELATION_IS_LOCAL(rel)))
			LockPage(rel, targetBlkno, ML_SHARE);
		buf = _ml_getbuf(rel, targetBlkno, ML_READ);
		npage = BufferGetPage(buf);

		Size modelsize = sizeof(MLModel);
		Size pgspc = PageGetFreeSpace(npage);

		 /* Check whether the model can fit on a mlree page at all.*/

		if(pgspc < modelsize){
			 /* get a new page.*/
			UnlockPage(rel, targetBlkno, ML_SHARE);
			targetBlkno = targetBlkno + 1;
			npage = _ml_newpage();

			LockPage(rel, targetBlkno, ML_SHARE);

			//page is up to MAXMETABLOCK,we alloc a new metapage to record meta info.
			if(targetBlkno % MAXMETABLOCK==0){
				//maybe here we don't need to mark the buffer dirty because it's dirty enough
				UnlockReleaseBuffer(metabuf);
				metabuf = _ml_metapinit(rel,targetBlkno++);
				metap = MLPageGetMeta(BufferGetPage(metabuf));
				MarkBufferDirty(metabuf);
			}

			metap->length++;
		}

		 /* Add the new item into the current page.*/


		MLPageAddItem(npage,(Item)mlmodel,modelsize);
		//length is wrong
		metap->pageInfo[metap->length-1].blockmax = mlmodel->maxvalue;
		metap->pageInfo[metap->length-1].blkno = mlmodel->blkno;
		_ml_writepage(rel, npage, targetBlkno,MAIN_FORKNUM);
		//INSTR_TIME_SET_CURRENT(after);
		//INSTR_TIME_SUBTRACT(after, before);
		//time1 += INSTR_TIME_GET_MILLISEC(after);
	}
	else{
		metablkno = _ml_findTargetPage(metap->pageInfo,blkno,metap->length);
		targetBlkno = i*MAXMETABLOCK+metablkno+1;
		/*
		 * Acquire shared split lock
		 */
		if ((!RELATION_IS_LOCAL(rel)))
			LockPage(rel, targetBlkno, ML_SHARE);
		/* Fetch the page*/
		buf = _ml_getbuf(rel, targetBlkno, ML_READ);
		page = BufferGetPage(buf);

		//INSTR_TIME_SET_CURRENT(before);
		//do binary search on targeted page
		MLModel* left = MLPageGetFirstModel(page);
		MLModel* right = MLPageGetLastModel(page);
		MLModel* mid = left + ((right - left)>>1);

		//find the proper model in the page
		while (left <= right)
		{
			mid = left + ((right - left)>>1);
			if (blkno <mid->blkno)
				right = mid-1;
			else if (blkno > mid->blkno)
				left = mid+1;
			else{
				break;
			}
		}

		targetModel = mid;
		/*
		Datum leftmax,rightmin;

		if(left == right){
			if(metablkno==0){
				if(i==0){
					leftmax = 0;
				}
				metabuf = _ml_getbuf(rel, (i-1)*MAXMETABLOCK, ML_READ);
				metap = MLPageGetMeta(BufferGetPage(metabuf));
				metablkno = metap->length;
			}
			if()

		}
		if(i == 0 && metablkno == 0){
			maxvalue = 0;
			if(left == right)
				minvalue = INT_MAX;
			else
				minvalue = (targetModel+1)->minvalue;
		}
		if(metap->isLastMeta && metablkno == metap->length-1){
			minvalue = INT_MAX;
			if(left == right)

			else
				minvalue = (targetModel+1)->minvalue;
		}
		if(metablkno==0){
			metabuf = _ml_getbuf(rel, (i-1)*MAXMETABLOCK, ML_READ);
			metap = MLPageGetMeta(BufferGetPage(metabuf));
			metablkno = metap->length;
		}
		else if(metablkno == metap->length-1){
			metabuf = _ml_getbuf(rel, (i+1)*MAXMETABLOCK, ML_READ);
			metap = MLPageGetMeta(BufferGetPage(metabuf));
			metablkno = -1;
		}
		if(targetModel == left){
			maxvalue = metap->pageInfo[metablkno-1].blockmax;
			minvalue = (targetModel+1)->minvalue;
		}
		else if(targetModel == right){
			maxvalue = (targetModel-1)->maxvalue;
			minvalue = metap->pageInfo[metablkno+1].blockmin;
		}
		else{
			maxvalue = (targetModel-1)->maxvalue;
			minvalue = (targetModel+1)->minvalue;
		}
		if(values[0] > maxvalue || values[0] < minvalue){

		}
*/
		//do some change to the model
		_ls_insert(values[0],pos,targetModel,heapRel,rel);

		metap->pageInfo[metap->length-1].blockmax = targetModel->maxvalue;

	}
	UnlockReleaseBuffer(buf);
	UnlockReleaseBuffer(metabuf);
	UnlockPage(rel, targetBlkno, ML_SHARE);
	INSTR_TIME_SET_CURRENT(after);
	INSTR_TIME_SUBTRACT(after, before);
	time1 += INSTR_TIME_GET_MILLISEC(after);
	PG_RETURN_BOOL(true);
}


Datum
mloptions(PG_FUNCTION_ARGS)
{

	PG_RETURN_NULL();
}

Datum
mlint4cmp(PG_FUNCTION_ARGS)
{
	int32		a = PG_GETARG_INT32(0);
	int32		b = PG_GETARG_INT32(1);

	if (a > b)
		PG_RETURN_INT32(1);
	else if (a == b)
		PG_RETURN_INT32(0);
	else
		PG_RETURN_INT32(-1);
}

