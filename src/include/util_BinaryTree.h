 /*@@
   @header    BinaryTree.h
   @date      Mon Oct  5 11:01:20 1998
   @author    Tom Goodale
   @desc 
   Prototypes and data definitions for binary tree routines.
   @enddesc 
   @version $Header$
 @@*/

#ifndef _BINARYTREE_H_
#define _BINARYTREE_H_ 1

#ifdef __cplusplus
extern "C" 
{
#endif

typedef struct T_TREE
{
  struct T_TREE *left;
  struct T_TREE *right;
  struct T_TREE *next;

  void *data;
} uBinTree;


uBinTree *Util_BinTreeStoreData(uBinTree *root, 
                                uBinTree *subtree, 
                                void *data, 
                                int (*compare)(const void *, const void *));

int Util_BinTreeTraverseInorder(uBinTree *root, 
                                int (*process)(void *, void *), 
                                void *info);

int Util_BinTreeTraversePreorder(uBinTree *root, 
                                 int (*process)(void *, void *), 
                                 void *info);

int Util_BinTreeTraversePostorder(uBinTree *root, 
                                  int (*process)(void *, void *), 
                                  void *info);

int Util_BinTreePrintNodes(uBinTree *root, 
                           int depth, 
                           void (*print_node)(void *, int));

uBinTree *Util_BinTreeFindNode(uBinTree *root, 
                               void *data, 
                               int (*compare)(const void *, const void *));

#ifdef __cplusplus
}
#endif

#endif
