//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Fri May 29 21:16:29 2026] 

#ifndef ROOT_TMVA_SOFIE_TRAINBKGDIAG
#define ROOT_TMVA_SOFIE_TRAINBKGDIAG

#include <algorithm>
#include <cmath>
#include <vector>
#include "TMVA/SOFIE_common.hxx"
#include <fstream>

namespace TMVA_SOFIE_TrainBkgDiag{
namespace BLAS{
	extern "C" void sgemv_(const char * trans, const int * m, const int * n, const float * alpha, const float * A,
	                       const int * lda, const float * X, const int * incx, const float * beta, const float * Y, const int * incy);
	extern "C" void sgemm_(const char * transa, const char * transb, const int * m, const int * n, const int * k,
	                       const float * alpha, const float * A, const int * lda, const float * B, const int * ldb,
	                       const float * beta, float * C, const int * ldc);
}//BLAS
struct Session {
std::vector<float> fTensor_dense15bias0 = std::vector<float>(1);
float * tensor_dense15bias0 = fTensor_dense15bias0.data();
std::vector<float> fTensor_dense15kernel0 = std::vector<float>(32);
float * tensor_dense15kernel0 = fTensor_dense15kernel0.data();
std::vector<float> fTensor_dense14bias0 = std::vector<float>(32);
float * tensor_dense14bias0 = fTensor_dense14bias0.data();
std::vector<float> fTensor_dense14kernel0 = std::vector<float>(2048);
float * tensor_dense14kernel0 = fTensor_dense14kernel0.data();
std::vector<float> fTensor_dense13bias0 = std::vector<float>(64);
float * tensor_dense13bias0 = fTensor_dense13bias0.data();
std::vector<float> fTensor_dense13kernel0 = std::vector<float>(4096);
float * tensor_dense13kernel0 = fTensor_dense13kernel0.data();
std::vector<float> fTensor_dense12bias0 = std::vector<float>(64);
float * tensor_dense12bias0 = fTensor_dense12bias0.data();
std::vector<float> fTensor_dense12kernel0 = std::vector<float>(448);
float * tensor_dense12kernel0 = fTensor_dense12kernel0.data();
std::vector<float> fTensor_dense14Relu0 = std::vector<float>(32);
float * tensor_dense14Relu0 = fTensor_dense14Relu0.data();
std::vector<float> fTensor_dense14Dense = std::vector<float>(32);
float * tensor_dense14Dense = fTensor_dense14Dense.data();
std::vector<float> fTensor_dense14bias0bcast = std::vector<float>(32);
float * tensor_dense14bias0bcast = fTensor_dense14bias0bcast.data();
std::vector<float> fTensor_dense13Relu0 = std::vector<float>(64);
float * tensor_dense13Relu0 = fTensor_dense13Relu0.data();
std::vector<float> fTensor_dense12Relu0 = std::vector<float>(64);
float * tensor_dense12Relu0 = fTensor_dense12Relu0.data();
std::vector<float> fTensor_dense15Dense = std::vector<float>(1);
float * tensor_dense15Dense = fTensor_dense15Dense.data();
std::vector<float> fTensor_dense15bias0bcast = std::vector<float>(1);
float * tensor_dense15bias0bcast = fTensor_dense15bias0bcast.data();
std::vector<float> fTensor_dense15Sigmoid0 = std::vector<float>(1);
float * tensor_dense15Sigmoid0 = fTensor_dense15Sigmoid0.data();
std::vector<float> fTensor_dense12bias0bcast = std::vector<float>(64);
float * tensor_dense12bias0bcast = fTensor_dense12bias0bcast.data();
std::vector<float> fTensor_dense13bias0bcast = std::vector<float>(64);
float * tensor_dense13bias0bcast = fTensor_dense13bias0bcast.data();
std::vector<float> fTensor_dense13Dense = std::vector<float>(64);
float * tensor_dense13Dense = fTensor_dense13Dense.data();
std::vector<float> fTensor_dense12Dense = std::vector<float>(64);
float * tensor_dense12Dense = fTensor_dense12Dense.data();


Session(std::string filename ="") {
   if (filename.empty()) filename = "TrainBkgDiag.dat";
   std::ifstream f;
   f.open(filename);
   if (!f.is_open()) {
      throw std::runtime_error("tmva-sofie failed to open file for input weights");
   }
   std::string tensor_name;
   size_t length;
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense15bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense15bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense15bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense15kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense15kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense15kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense14bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense14bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense14bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense14kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense14kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 2048) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 2048 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense14kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense13bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense13bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense13bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense13kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense13kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense13kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense12bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense12bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense12bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense12kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense12kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 448) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 448 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense12kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense12bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense12bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense13bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense13bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense14bias0,{ 32 }, { 1 , 32 });
      std::copy(data, data + 32, tensor_dense14bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense15bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense15bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input4){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 1;
   int op_0_n = 64;
   int op_0_k = 7;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 7;
   int op_0_ldb = 64;
   std::copy(tensor_dense12bias0bcast, tensor_dense12bias0bcast + 64, tensor_dense12Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense12kernel0, &op_0_ldb, tensor_input4, &op_0_lda, &op_0_beta, tensor_dense12Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense12Relu0[id] = ((tensor_dense12Dense[id] > 0 )? tensor_dense12Dense[id] : 0);
   }

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 1;
   int op_2_n = 64;
   int op_2_k = 64;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 64;
   int op_2_ldb = 64;
   std::copy(tensor_dense13bias0bcast, tensor_dense13bias0bcast + 64, tensor_dense13Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense13kernel0, &op_2_ldb, tensor_dense12Relu0, &op_2_lda, &op_2_beta, tensor_dense13Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense13Relu0[id] = ((tensor_dense13Dense[id] > 0 )? tensor_dense13Dense[id] : 0);
   }

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 1;
   int op_4_n = 32;
   int op_4_k = 64;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 64;
   int op_4_ldb = 32;
   std::copy(tensor_dense14bias0bcast, tensor_dense14bias0bcast + 32, tensor_dense14Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense14kernel0, &op_4_ldb, tensor_dense13Relu0, &op_4_lda, &op_4_beta, tensor_dense14Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 32 ; id++){
      tensor_dense14Relu0[id] = ((tensor_dense14Dense[id] > 0 )? tensor_dense14Dense[id] : 0);
   }

//--------- Gemm
   char op_6_transA = 'n';
   char op_6_transB = 'n';
   int op_6_m = 1;
   int op_6_n = 1;
   int op_6_k = 32;
   float op_6_alpha = 1;
   float op_6_beta = 1;
   int op_6_lda = 32;
   int op_6_ldb = 1;
   std::copy(tensor_dense15bias0bcast, tensor_dense15bias0bcast + 1, tensor_dense15Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense15kernel0, &op_6_ldb, tensor_dense14Relu0, &op_6_lda, &op_6_beta, tensor_dense15Dense, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense15Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense15Dense[id]));
	}
   std::vector<float> ret (tensor_dense15Sigmoid0, tensor_dense15Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
