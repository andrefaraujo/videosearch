// -*- c++ -*- 


/**************************************************************
 * Useful stuff to include in other SWIG interfaces:
 * - macros to wrap slow functions that need to release the GIL 
 * - FILE * typemap to a python file object
 */

%{

/* To call a peek function from a slow function that released the GIL:
 * GIL must be acquired first, then callback, then released again
 */
typedef struct {
  PyObject *callable;
  PyThreadState *_save;
} callback_t;

static void call_peek_fun_thread(void *arg,double frac) {
  callback_t *callback=(callback_t*)arg;
  PyEval_RestoreThread(callback->_save);  
  {
    PyObject *obj=PyObject_CallFunction(callback->callable,"d",frac);
    if(!obj) {
      fprintf(stderr,"!!!exception in peek callback\n");
      PyErr_Print();
    } else
      Py_DECREF(obj);   
  }
  callback->_save=PyEval_SaveThread();
}

/* cancels renaming of callback_t *cb */
#define cb1 cb
#define cb2 cb
#define cb3 cb
#define cb4 cb
#define cb5 cb
#define cb6 cb
#define cb7 cb
#define cb8 cb
#define cb9 cb
#define cb10 cb
#define cb11 cb
#define cb12 cb
#define cb13 cb


%}


// function releases GIL while doing slow action
%define PARALLEL_FUNCTION(funcname)
%exception funcname {
  Py_BEGIN_ALLOW_THREADS
  $action
  Py_END_ALLOW_THREADS
}
%enddef

/*
 * same as before, but must be locked from time to time to 
 * call peek callback, hence the callback_t structure that holds both the 
 * function to call and the GIL. 
 */

%define PARALLEL_FUNCTION_PEEK(funcname)
%exception funcname {
  cb._save = PyEval_SaveThread();
  $action
  PyEval_RestoreThread(cb._save);
}
%enddef

%typemap(in, numinputs=1) 
    (void (*peek_fun)(void *arg,double frac),void *peek_arg) 
    (callback_t cb)
{
  if(Py_None==$input) {
    $1=NULL;
    $2=NULL; 
  } else {
    cb.callable=$input;
    $1=call_peek_fun_thread;
    $2=&cb;
  }
}



// map Python file objects (obtained with open()) to FILE*

%typemap(in) FILE* {
  if(!PyFile_Check($input)) {
    PyErr_SetString(PyExc_TypeError,"need regular file object");
    return NULL;
  }
  $1=((PyFileObject*)$input)->f_fp;
}


%include <carrays.i>

/* 
   a=IntArray(4) constructs an array of 4 ints  
   a.cast() returns an int* usable as a function argument
   if x is a int* in C,  
   b=IntArray.frompointer(x) makes a[i] valid to access x[i]
   b=IntArray.acquirepointer(x) will also free(x) when b is deleted
   b.plus(2) returns a pointer to x+2
*/
%define array_class_plus(type,name)
%array_class(type,name)
%extend name {
  type *plus(long i) {
    return $self+i;
  }
  void clear(long n) {
    memset($self,0,sizeof(type)*n);
  }
  void copyfrom(type *src,long dest_ofs,long n) {
    memmove($self+dest_ofs,src,sizeof(type)*n);
  } 
  PyObject *tostring(long n) {
    PyObject *str = PyString_FromStringAndSize(NULL, n*sizeof(type));
    memcpy(PyString_AsString(str),$self,sizeof(type)*n);
    return str;
  }
  void fromstring(PyObject*obj) {
    if(!PyString_Check(obj)) {
      fprintf(stderr,"need string!");
      return;
    }
    memcpy($self,PyString_AsString(obj),PyString_Size(obj));
  }
};
%pythoncode %{
name.aptr=name.acquirepointer=staticmethod(lambda ptr: _frompointer_and_acquire(name,ptr))
%}
%enddef

%pythoncode %{
def _frompointer_and_acquire(aclass,ptr):
  r=aclass.frompointer(ptr)
  if r: r.this.acquire()
  return r
%}
