// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME G__LUXESimEventDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// Header files passed as explicit arguments
#include "./include/ActsLUXEPipeline/detail/ROOTDictDefs.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *vectorlETVector3gR_Dictionary();
   static void vectorlETVector3gR_TClassManip(TClass*);
   static void *new_vectorlETVector3gR(void *p = nullptr);
   static void *newArray_vectorlETVector3gR(Long_t size, void *p);
   static void delete_vectorlETVector3gR(void *p);
   static void deleteArray_vectorlETVector3gR(void *p);
   static void destruct_vectorlETVector3gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector3>*)
   {
      vector<TVector3> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector3>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector3>", -2, "vector", 383,
                  typeid(vector<TVector3>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETVector3gR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TVector3>) );
      instance.SetNew(&new_vectorlETVector3gR);
      instance.SetNewArray(&newArray_vectorlETVector3gR);
      instance.SetDelete(&delete_vectorlETVector3gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector3gR);
      instance.SetDestructor(&destruct_vectorlETVector3gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector3> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TVector3>","std::__1::vector<TVector3, std::__1::allocator<TVector3>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TVector3>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector3gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TVector3>*>(nullptr))->GetClass();
      vectorlETVector3gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector3gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector3gR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TVector3> : new vector<TVector3>;
   }
   static void *newArray_vectorlETVector3gR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TVector3>[nElements] : new vector<TVector3>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector3gR(void *p) {
      delete (static_cast<vector<TVector3>*>(p));
   }
   static void deleteArray_vectorlETVector3gR(void *p) {
      delete [] (static_cast<vector<TVector3>*>(p));
   }
   static void destruct_vectorlETVector3gR(void *p) {
      typedef vector<TVector3> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TVector3>

namespace ROOT {
   static TClass *vectorlETLorentzVectorgR_Dictionary();
   static void vectorlETLorentzVectorgR_TClassManip(TClass*);
   static void *new_vectorlETLorentzVectorgR(void *p = nullptr);
   static void *newArray_vectorlETLorentzVectorgR(Long_t size, void *p);
   static void delete_vectorlETLorentzVectorgR(void *p);
   static void deleteArray_vectorlETLorentzVectorgR(void *p);
   static void destruct_vectorlETLorentzVectorgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TLorentzVector>*)
   {
      vector<TLorentzVector> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TLorentzVector>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TLorentzVector>", -2, "vector", 383,
                  typeid(vector<TLorentzVector>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETLorentzVectorgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<TLorentzVector>) );
      instance.SetNew(&new_vectorlETLorentzVectorgR);
      instance.SetNewArray(&newArray_vectorlETLorentzVectorgR);
      instance.SetDelete(&delete_vectorlETLorentzVectorgR);
      instance.SetDeleteArray(&deleteArray_vectorlETLorentzVectorgR);
      instance.SetDestructor(&destruct_vectorlETLorentzVectorgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TLorentzVector> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<TLorentzVector>","std::__1::vector<TLorentzVector, std::__1::allocator<TLorentzVector>>"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<TLorentzVector>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETLorentzVectorgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<TLorentzVector>*>(nullptr))->GetClass();
      vectorlETLorentzVectorgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETLorentzVectorgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETLorentzVectorgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TLorentzVector> : new vector<TLorentzVector>;
   }
   static void *newArray_vectorlETLorentzVectorgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<TLorentzVector>[nElements] : new vector<TLorentzVector>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETLorentzVectorgR(void *p) {
      delete (static_cast<vector<TLorentzVector>*>(p));
   }
   static void deleteArray_vectorlETLorentzVectorgR(void *p) {
      delete [] (static_cast<vector<TLorentzVector>*>(p));
   }
   static void destruct_vectorlETLorentzVectorgR(void *p) {
      typedef vector<TLorentzVector> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<TLorentzVector>

namespace {
  void TriggerDictionaryInitialization_libLUXESimEventDict_Impl() {
    static const char* headers[] = {
"./include/ActsLUXEPipeline/detail/ROOTDictDefs.h",
nullptr
    };
    static const char* includePaths[] = {
"/Users/alonlevi/CLionProjects/LUXEPipeline/Readout/include/ActsLUXEPipeline/detail",
"/Users/alonlevi/CLionProjects/LUXEPipeline/Readout",
"/Users/alonlevi/CLionProjects/LUXEPipeline/Readout/include/ActsLUXEPipeline/detail",
"/Users/alonlevi/Programs/root_install/include",
"/Users/alonlevi/Programs/root_install/include/",
"/Users/alonlevi/CLionProjects/LUXEPipeline/build/Readout/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libLUXESimEventDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$TVector3.h")))  __attribute__((annotate("$clingAutoload$./include/ActsLUXEPipeline/detail/ROOTDictDefs.h")))  TVector3;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$__fwd/string.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
class __attribute__((annotate("$clingAutoload$TLorentzVector.h")))  __attribute__((annotate("$clingAutoload$./include/ActsLUXEPipeline/detail/ROOTDictDefs.h")))  TLorentzVector;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libLUXESimEventDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./include/ActsLUXEPipeline/detail/ROOTDictDefs.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libLUXESimEventDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libLUXESimEventDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libLUXESimEventDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libLUXESimEventDict() {
  TriggerDictionaryInitialization_libLUXESimEventDict_Impl();
}
