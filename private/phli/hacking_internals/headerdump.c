/* headerdump.c: MEX file to show the internals of an mxArray
   header.

   WARNING! WARNING! WARNING!  This program blatantly abuses
   the Matlab API!  It pokes around in private data structures
   without remorse!  The structures may change from version to
   version (or even platform to platform for all we know).

   That said, this program only looks at memory.  So even if you
   get a Matlab segfault, nothing should be corrupted.  Please
   let me know, though, if you do segfault, with instructions on
   how to reproduce it.

   The author and his employer disclaim all liability for any damage
   caused by this program.  For experimental use only!

   Author: Peter Boettcher <boettcher@ll.mit.edu>
   Last modified: <Mon Jul 14 11:05:38 2003 by pwb> 

   Support for R13 added by Gareth Jones.  Thanks!
*/

#include "mex.h"
#include "mxinternals.h"

bool mxIsSharedArray(const mxArray *array_ptr);
void mxUnshareArray(mxArray *array_ptr);

/* Print one-liner describing mxArray, including any crosslinks */
void briefdump(mxArray *in)
{
  int *tmp;
  int i;

  if (in->name && *in->name)
    printf("[%.32s]", in->name);
  printf("(%8s)  ", mxGetClassName(in));
  printf("Address: %p", in);
  if(in->crosslink)
    printf("  <linked to %p>", in->crosslink);

  if(in->refcount)
    printf("  <refcount %d>", in->refcount);

  if(in->number_of_dims == 2) {
    printf("   [%i %i]", in->rowdim, in->coldim);
    if(in->rowdim*in->coldim == 1 && mxGetClassID(in)==mxDOUBLE_CLASS)
      printf(" scalar value: %f", *((double*)(in->data.number_array.pdata)));
    printf("\n");
  } else {
    printf("   [");
    tmp = (int *)(in->rowdim);
    for(i=0; i<in->number_of_dims; i++) {
      printf("%i ", tmp[i]);
    }
    printf("]");
  }
}

/* Print detailed description of mxArray */
void dumpMxArray(const mxArray *in)
{
  int *tmp;
  int i;
  int numel;
  void *tmpp;

  printf("Name: %.32s\n", (in->name && *in->name) ? in->name : "<temp>");
  printf("Address: %p", in);

  /* Crosslink means two or more variables point to the same
     data.  The link allows Matlab to copy the array and update
     the affected variables if someone wants to modify an element */
  if(in->crosslink)
    printf("   <Crosslinked to %p: %.32s>\n", in->crosslink,
	   in->crosslink->name);
  else
    printf("\n");

#if MATLAB_RELEASE > 11
  if(in->refcount)
    printf("   Refcount: %d\n", in->refcount);
#endif

  /* This field has a unique value for each class, but is not
     the same as the class ID */
  printf("Related to classID?  %i (true: %i)\n",
	 in->class_id, (int)mxGetClassID(in));
  printf("%s\n", mxGetClassName(in));
  printf("Variable type: ");
  switch(in->vartype) {
  case MXVARNORMAL:
    printf("Normal\n");
    break;
  case MXVARPERSIST:
    printf("Persistent\n");
    break;
  case MXVARGLOBAL:
    printf("Global\n");
    break;
  case MXVARSUBEL:
    printf("Subelement of a cell or struct\n");
    break;
  case MXVARTEMP:
    printf("Temporary\n");
    break;
  default:
    printf("Unknown variable type: %i  (Please email boettcher@ll.mit.edu)\n",
	   in->vartype);
  }

#if MATLAB_RELEASE < 13
  printf("Data Flags: Numeric? %i Empty %i  Logical %i  DblScalar %i  (other: %x)\n",
	 (in->dataflags&MXNUMERICMASK)!=0,
	 (in->dataflags&MXEMPTYMASK)!=0,
	 (in->dataflags&MXLOGICALMASK)!=0,
	 (in->dataflags&MXSCALARMASK)!=0,
	 (in->dataflags&0x00fffdf8));
  if(in->dataflags & 0xff000000)
    printf("User Data: %x\n", (in->dataflags&0xff000000)>>24);

#else
  printf("Data Flags: DblScalar %i  private_data_flag %i\n",
	 in->flags.scalar_flag!=0,
	 in->flags.private_data_flag!=0);
#endif

  printf("\nDimensions (%i): ", in->number_of_dims);
  if(in->number_of_dims == 2) {
    printf("[%i %i]\n", in->rowdim, in->coldim);
    numel = in->rowdim * in->coldim;
  } else { /* multidimensional */
    numel = 1;
    printf("[");
    tmp = (int *)(in->rowdim);
    for(i=0; i<in->number_of_dims; i++) {
      numel *= tmp[i];
      printf("%i ", tmp[i]);
    }
    printf("]\n");
  }
  
  printf("Real data: %p", in->data.number_array.pdata);
#if MATLAB_RELEASE < 13
  if(in->dataflags & MXSCALARMASK)
#else /*MATLAB_RELEASE == 13*/
  if(in->flags.scalar_flag)
#endif
     printf(" [%g]", *((double *)in->data.number_array.pdata));
  printf("\n");
  
  if(in->data.number_array.pimag_data) { /* complex */
    printf("Imag data: %p", in->data.number_array.pimag_data);
#if MATLAB_RELEASE < 13
    if(in->dataflags & MXSCALARMASK)
#else
    if(in->flags.scalar_flag)
#endif
      printf(" [%g]", *((double *)in->data.number_array.pimag_data));
    printf("\n");
  }

#if MATLAB_RELEASE > 11
//  if(in->nelements_allocated) { /* sparse */
//    printf("Sparse matrix: Nelements: %i\n", in->nelements_allocated);
//    printf("Column ptr: %p\nRow ptr: %p\n",
//	   in->data.number_array.irptr, in->data.number_array.jcptr);
//  }
#else
  if(in->class_id == 3) {
    printf("Sparse matrix: Nelements: %i\n", in->data.number_array.nelements);
    printf("Column ptr: %p\nRow ptr: %p\n",
	   in->data.number_array.irptr, in->data.number_array.jcptr);
  }
  /* else if (in->data.number_array.nelements ||	in->data.number_array.irptr 
	     ||in->data.number_array.jcptr)
    printf("Not sparse BUT! : %p %p %p\n", in->data.number_array.nelements,
	   in->data.number_array.irptr, in->data.number_array.jcptr);
  */
  else if(in->class_id == 18) {
    printf("Object:  Name: %s\n", in->data.object_array.name);
    printf("         Checksum??: 0x%x\n", in->data.object_array.checksum);
    printf("         Nelements??: 0x%x\n", in->data.object_array.nelements);
  }

  else if(in->class_id == 17) {
    for(i=0; i<10; i++) {
      printf("0x%x ", *(((unsigned long *)in->data.number_array.pdata) + i));
    }
    printf("\n");
    tmpp = *((char **)in->data.number_array.pdata + 1);
    for(i=0; i<10; i++) {
      printf("0x%x ", *(((void **)tmpp) + i));
    }
    printf("\n");
    tmpp = *((void**)tmpp);
    printf("%s\n", (char *)tmpp);

    tmpp = *((char **)in->data.number_array.pdata);
    for(i=0; i<10; i++) {
      printf("0x%x ", *(((void **)tmpp) + i));
    }
    printf("\n");
  }
#endif
  
#if 0
  if(in->data.number_array.reserved) /* what's this? */ {
    printf("Unknown field: %x (Please email boettcher@ll.mit.edu)\n",
	   in->data.number_array.reserved);
  }
#endif

  if(in->class_id == 6 || in->class_id == 18) /* struct */ {
    printf("Number of struct fields: %i\n", in->data.number_array.nfields);
    numel *= in->data.struct_array.nfields;
    for(i=0; i<in->data.struct_array.nfields; i++)
      printf(" %.32s\n", in->data.struct_array.field_names + i*32);
    printf("\n");
  }
  printf("\n");

  /* Structs are stored like cells, using mxArray pointers in the real
     data spot.  They are stored in "field major" order, meaning the
     pointers to the values of the fields of the first struct element
     are stored in order, then the columns of the struct array, then
     the rows. */
  /* Print a one-liner for each element of the cell or struct, up
     to a maximum of 30 */
  if(in->class_id == 6 || in->class_id == 5 || in->class_id == 18) {
    for(i=0; i<((numel < 30) ? numel : 30); i++)
      if(((mxArray **)in->data.number_array.pdata)[i])
	briefdump(((mxArray **)in->data.number_array.pdata)[i]);
      else
	printf("(nil)\n");
  }
}  


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if(nrhs < 1) 
    mexErrMsgTxt("One input required.");  

  dumpMxArray(prhs[0]);
}
