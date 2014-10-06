/*
 * Reverse engineered header for Matlab's internal mxArray structure.
 * 
 * Author: Peter Boettcher <boettcher@ll.mit.edu>
 *
 * Support for R13 added by Gareth Jones.  Thanks!
 */

#ifndef MXINTERNALS_H
#define MXINTERNALS_H 1

#define MATLAB_RELEASE 13


#define MXVARNORMAL (0)
#define MXVARPERSIST (1)
#define MXVARGLOBAL (2)
#define MXVARSUBEL (3)
#define MXVARTEMP (4)

#if 0
#define MXKERNELBITS 0xf000


#define MXCOMMALISTMASK 0x800  /* mxSetCommaList */
#define MXDATAPRIVATE 0x100 /* mxSetDataPrivate */
#endif
#define MXNUMERICMASK 0x200
#define MXCOLONMASK 0x08
#define MXEMPTYMASK 0x04
#define MXLOGICALMASK 0x02
#define MXSCALARMASK 0x01

/* Blatant disregard for MATLAB API.  This allows you
   to poke around an mxArray structure directly.  Use
   with extreme caution!  May not be portable across
   platforms or versions, may cause Matlab segfaults,
   may cause your computer to explode and eat all your data... */
#if MATLAB_RELEASE == 11
struct mxArray_tag {
  char     name[mxMAXNAM];
  int class_id;
  int vartype;
  mxArray    *crosslink;
  int      number_of_dims;
  int      nelements_allocated;
  int  dataflags;
  int  rowdim;
  int  coldim;
  union {
    struct {
      void  *pdata;
      void  *pimag_data;
      void  *irptr;
      void  *jcptr;
      int    reserved;
      int   nfields;
    }   number_array;
  }   data;
};
#elif MATLAB_RELEASE == 12
struct mxArray_tag {
  char     name[mxMAXNAM];
  int class_id;
  int vartype;
  mxArray    *crosslink;
  int      number_of_dims;
  int      refcount;
  int  dataflags;
  int  rowdim;
  int  coldim;
  union {
    struct {
      void  *pdata;
      void  *pimag_data;
      void  *irptr;
      void  *jcptr;
      int   nelements;
      int   nfields;
    }   number_array;
    struct {
      mxArray **pdata;
      char  *field_names;
      void  *dummy1;
      void  *dummy2;
      int   dummy3;
      int   nfields;
    }   struct_array;
    struct {
      void *pdata;  /*mxGetInfo*/
      void *pimag_data;
      void *irptr;
      void *jcptr;
      int  nelements;
      int  reserved;
    }  sparse_array;
    struct {
      void *pdata;  /*mxGetInfo*/
      char *field_names;
      char *name;
      int checksum;
      int  nelements;
      int  reserved;
    }  object_array;
  }   data;
};
#elif MATLAB_RELEASE == 13
struct mxArray_tag {
  const char *name;
  int class_id;
  int vartype;
  mxArray    *crosslink;
  int      number_of_dims;
  int      refcount;
  struct {
    unsigned int    scalar_flag : 1;
    unsigned int	flag1 : 1; 
    unsigned int    flag2 : 1;
    unsigned int    flag3 : 1;
    unsigned int    flag4 : 1;
    unsigned int    flag5 : 1;
    unsigned int    flag6 : 1;
    unsigned int    flag7 : 1;
    unsigned int    private_data_flag : 1;
    unsigned int    flag8 : 1;
    unsigned int    flag9 : 1;
    unsigned int    flag10 : 1;
    unsigned int    flag11 : 4;
    unsigned int    flag12 : 8;
    unsigned int    flag13 : 8;
  }   flags;
  int  rowdim;
  int  coldim;
  union {
    struct {
      void  *pdata;
      void  *pimag_data;
      void *irptr;
      void  *jcptr;
      int   nelements;
      int   nfields;
    }   number_array;
    struct {
      mxArray **pdata;
      char  *field_names;
      void  *dummy1;
      void  *dummy2;
      int   dummy3;
      int   nfields;
    }   struct_array;
    struct {
      void *pdata;  /*mxGetInfo*/
      char *field_names;
      char *name;
      int checksum;
      int  nelements;
      int  reserved;
    }  object_array;
  }   data;
};
#elif MATLAB_RELEASE == R2010a
// Copied from extern/include/matrix.h
struct mxArray_tag {
    void    *reserved;
    int      reserved1[2];
    void    *reserved2;
    size_t  number_of_dims;
    unsigned int reserved3;
    struct {
        unsigned int    flag0 : 1;
        unsigned int    flag1 : 1;
        unsigned int    flag2 : 1;
        unsigned int    flag3 : 1;
        unsigned int    flag4 : 1;
        unsigned int    flag5 : 1;
        unsigned int    flag6 : 1;
        unsigned int    flag7 : 1;
        unsigned int    flag7a: 1;
        unsigned int    flag8 : 1;
        unsigned int    flag9 : 1;
        unsigned int    flag10 : 1;
        unsigned int    flag11 : 4;
        unsigned int    flag12 : 8;
        unsigned int    flag13 : 8;
    }   flags;
    size_t reserved4[2];
    union {
        struct {
            void  *pdata;
            void  *pimag_data;
            void  *reserved5;
            size_t reserved6[3];
        }   number_array;
    }   data;
};

#else
#error Unknown matlab release MATLAB_RELEASE
#endif

#endif
