// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME LoadInputDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "LoadInput.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_Hit(void *p = 0);
   static void *newArray_Hit(Long_t size, void *p);
   static void delete_Hit(void *p);
   static void deleteArray_Hit(void *p);
   static void destruct_Hit(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Hit*)
   {
      ::Hit *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Hit >(0);
      static ::ROOT::TGenericClassInfo 
         instance("Hit", ::Hit::Class_Version(), "LoadInput.h", 10,
                  typeid(::Hit), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Hit::Dictionary, isa_proxy, 4,
                  sizeof(::Hit) );
      instance.SetNew(&new_Hit);
      instance.SetNewArray(&newArray_Hit);
      instance.SetDelete(&delete_Hit);
      instance.SetDeleteArray(&deleteArray_Hit);
      instance.SetDestructor(&destruct_Hit);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Hit*)
   {
      return GenerateInitInstanceLocal((::Hit*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::Hit*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SRawEvent(void *p = 0);
   static void *newArray_SRawEvent(Long_t size, void *p);
   static void delete_SRawEvent(void *p);
   static void deleteArray_SRawEvent(void *p);
   static void destruct_SRawEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SRawEvent*)
   {
      ::SRawEvent *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SRawEvent >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SRawEvent", ::SRawEvent::Class_Version(), "LoadInput.h", 25,
                  typeid(::SRawEvent), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SRawEvent::Dictionary, isa_proxy, 4,
                  sizeof(::SRawEvent) );
      instance.SetNew(&new_SRawEvent);
      instance.SetNewArray(&newArray_SRawEvent);
      instance.SetDelete(&delete_SRawEvent);
      instance.SetDeleteArray(&deleteArray_SRawEvent);
      instance.SetDestructor(&destruct_SRawEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SRawEvent*)
   {
      return GenerateInitInstanceLocal((::SRawEvent*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SRawEvent*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Hit::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Hit::Class_Name()
{
   return "Hit";
}

//______________________________________________________________________________
const char *Hit::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Hit::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Hit::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Hit::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Hit*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SRawEvent::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SRawEvent::Class_Name()
{
   return "SRawEvent";
}

//______________________________________________________________________________
const char *SRawEvent::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SRawEvent*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SRawEvent::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SRawEvent*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SRawEvent::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SRawEvent*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SRawEvent::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SRawEvent*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void Hit::Streamer(TBuffer &R__b)
{
   // Stream an object of class Hit.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Hit::Class(),this);
   } else {
      R__b.WriteClassBuffer(Hit::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Hit(void *p) {
      return  p ? new(p) ::Hit : new ::Hit;
   }
   static void *newArray_Hit(Long_t nElements, void *p) {
      return p ? new(p) ::Hit[nElements] : new ::Hit[nElements];
   }
   // Wrapper around operator delete
   static void delete_Hit(void *p) {
      delete ((::Hit*)p);
   }
   static void deleteArray_Hit(void *p) {
      delete [] ((::Hit*)p);
   }
   static void destruct_Hit(void *p) {
      typedef ::Hit current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::Hit

//______________________________________________________________________________
void SRawEvent::Streamer(TBuffer &R__b)
{
   // Stream an object of class SRawEvent.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SRawEvent::Class(),this);
   } else {
      R__b.WriteClassBuffer(SRawEvent::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SRawEvent(void *p) {
      return  p ? new(p) ::SRawEvent : new ::SRawEvent;
   }
   static void *newArray_SRawEvent(Long_t nElements, void *p) {
      return p ? new(p) ::SRawEvent[nElements] : new ::SRawEvent[nElements];
   }
   // Wrapper around operator delete
   static void delete_SRawEvent(void *p) {
      delete ((::SRawEvent*)p);
   }
   static void deleteArray_SRawEvent(void *p) {
      delete [] ((::SRawEvent*)p);
   }
   static void destruct_SRawEvent(void *p) {
      typedef ::SRawEvent current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SRawEvent

namespace ROOT {
   static TClass *vectorlEHitgR_Dictionary();
   static void vectorlEHitgR_TClassManip(TClass*);
   static void *new_vectorlEHitgR(void *p = 0);
   static void *newArray_vectorlEHitgR(Long_t size, void *p);
   static void delete_vectorlEHitgR(void *p);
   static void deleteArray_vectorlEHitgR(void *p);
   static void destruct_vectorlEHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<Hit>*)
   {
      vector<Hit> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<Hit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<Hit>", -2, "vector", 386,
                  typeid(vector<Hit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEHitgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<Hit>) );
      instance.SetNew(&new_vectorlEHitgR);
      instance.SetNewArray(&newArray_vectorlEHitgR);
      instance.SetDelete(&delete_vectorlEHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEHitgR);
      instance.SetDestructor(&destruct_vectorlEHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<Hit> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<Hit>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<Hit>*)0x0)->GetClass();
      vectorlEHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Hit> : new vector<Hit>;
   }
   static void *newArray_vectorlEHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<Hit>[nElements] : new vector<Hit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEHitgR(void *p) {
      delete ((vector<Hit>*)p);
   }
   static void deleteArray_vectorlEHitgR(void *p) {
      delete [] ((vector<Hit>*)p);
   }
   static void destruct_vectorlEHitgR(void *p) {
      typedef vector<Hit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<Hit>

namespace {
  void TriggerDictionaryInitialization_LoadInputDict_Impl() {
    static const char* headers[] = {
"LoadInput.h",
0
    };
    static const char* includePaths[] = {
"/home/cayuso/Products/test_build_6.14/include",
"/home/cayuso/Documents/e1039-GPU-master_current_copy/OnlineReconstruction/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "LoadInputDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$LoadInput.h")))  Hit;
class __attribute__((annotate("$clingAutoload$LoadInput.h")))  SRawEvent;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "LoadInputDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "LoadInput.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"Hit", payloadCode, "@",
"SRawEvent", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("LoadInputDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_LoadInputDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_LoadInputDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_LoadInputDict() {
  TriggerDictionaryInitialization_LoadInputDict_Impl();
}
