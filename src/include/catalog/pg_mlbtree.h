/*
 * pg_mlbtree.h
 *
 *  Created on: Oct 22, 2018
 *      Author: gpadmin
 */

#ifndef PG_MLBTREE_H
#define PG_MLBTREE_H

/*-------------------------------------------------------------------------
 */


#include "catalog/genbki.h"

/* ----------------
 * ----------------
 */
#define MLBTREERelationId	12600

CATALOG(pg_mlbtree,12600)
{
	Oid mltree;
	Oid btree;
} FormData_pg__mlbtree;

/* ----------------
 *		Form_pg_am corresponds to a pointer to a tuple with
 *		the format of pg_am relation.
 * ----------------
 */
typedef FormData_pg__mlbtree *Form_pg_mlbtree;

/* ----------------
 *		compiler constants for pg_am
 * ----------------
 */
#define Natts_pg_mlbtree					2
#define Anum_pg_mlbtree_ml				1
#define Anum_pg_mlbtree_btree				2



#endif   /* PG_MLBTREE_H */

