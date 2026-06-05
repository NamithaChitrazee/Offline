//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Fri Jun  5 18:48:44 2026] 

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
std::vector<float> fTensor_dense19bias0 = std::vector<float>(1);
float * tensor_dense19bias0 = fTensor_dense19bias0.data();
std::vector<float> fTensor_dense18bias0 = std::vector<float>(32);
float * tensor_dense18bias0 = fTensor_dense18bias0.data();
std::vector<float> fTensor_dense16bias0 = std::vector<float>(64);
float * tensor_dense16bias0 = fTensor_dense16bias0.data();
std::vector<float> fTensor_dense18kernel0 = std::vector<float>(2048);
float * tensor_dense18kernel0 = fTensor_dense18kernel0.data();
std::vector<float> fTensor_dense19kernel0 = std::vector<float>(32);
float * tensor_dense19kernel0 = fTensor_dense19kernel0.data();
std::vector<float> fTensor_dense17bias0 = std::vector<float>(64);
float * tensor_dense17bias0 = fTensor_dense17bias0.data();
std::vector<float> fTensor_dense17kernel0 = std::vector<float>(4096);
float * tensor_dense17kernel0 = fTensor_dense17kernel0.data();
std::vector<float> fTensor_dense16kernel0 = std::vector<float>(640);
float * tensor_dense16kernel0 = fTensor_dense16kernel0.data();
std::vector<float> fTensor_dense17Dense = std::vector<float>(64);
float * tensor_dense17Dense = fTensor_dense17Dense.data();
std::vector<float> fTensor_dense18Dense = std::vector<float>(32);
float * tensor_dense18Dense = fTensor_dense18Dense.data();
std::vector<float> fTensor_dense17Relu0 = std::vector<float>(64);
float * tensor_dense17Relu0 = fTensor_dense17Relu0.data();
std::vector<float> fTensor_dense18Relu0 = std::vector<float>(32);
float * tensor_dense18Relu0 = fTensor_dense18Relu0.data();
std::vector<float> fTensor_dense17bias0bcast = std::vector<float>(64);
float * tensor_dense17bias0bcast = fTensor_dense17bias0bcast.data();
std::vector<float> fTensor_dense18bias0bcast = std::vector<float>(32);
float * tensor_dense18bias0bcast = fTensor_dense18bias0bcast.data();
std::vector<float> fTensor_dense16Dense = std::vector<float>(64);
float * tensor_dense16Dense = fTensor_dense16Dense.data();
std::vector<float> fTensor_dense19Dense = std::vector<float>(1);
float * tensor_dense19Dense = fTensor_dense19Dense.data();
std::vector<float> fTensor_dense16Relu0 = std::vector<float>(64);
float * tensor_dense16Relu0 = fTensor_dense16Relu0.data();
std::vector<float> fTensor_dense19bias0bcast = std::vector<float>(1);
float * tensor_dense19bias0bcast = fTensor_dense19bias0bcast.data();
std::vector<float> fTensor_dense19Sigmoid0 = std::vector<float>(1);
float * tensor_dense19Sigmoid0 = fTensor_dense19Sigmoid0.data();
std::vector<float> fTensor_dense16bias0bcast = std::vector<float>(64);
float * tensor_dense16bias0bcast = fTensor_dense16bias0bcast.data();


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
   if (tensor_name != "tensor_dense19bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense19bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense19bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense18bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense18bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense18bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense16bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense16bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense16bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense18kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense18kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 2048) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 2048 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense18kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense19kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense19kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense19kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense17bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense17bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense17bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense17kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense17kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense17kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense16kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense16kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 640) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 640 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense16kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense16bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense16bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense17bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense17bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense18bias0,{ 32 }, { 1 , 32 });
      std::copy(data, data + 32, tensor_dense18bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense19bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense19bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input5){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 1;
   int op_0_n = 64;
   int op_0_k = 10;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 10;
   int op_0_ldb = 64;
   std::copy(tensor_dense16bias0bcast, tensor_dense16bias0bcast + 64, tensor_dense16Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense16kernel0, &op_0_ldb, tensor_input5, &op_0_lda, &op_0_beta, tensor_dense16Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense16Relu0[id] = ((tensor_dense16Dense[id] > 0 )? tensor_dense16Dense[id] : 0);
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
   std::copy(tensor_dense17bias0bcast, tensor_dense17bias0bcast + 64, tensor_dense17Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense17kernel0, &op_2_ldb, tensor_dense16Relu0, &op_2_lda, &op_2_beta, tensor_dense17Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense17Relu0[id] = ((tensor_dense17Dense[id] > 0 )? tensor_dense17Dense[id] : 0);
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
   std::copy(tensor_dense18bias0bcast, tensor_dense18bias0bcast + 32, tensor_dense18Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense18kernel0, &op_4_ldb, tensor_dense17Relu0, &op_4_lda, &op_4_beta, tensor_dense18Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 32 ; id++){
      tensor_dense18Relu0[id] = ((tensor_dense18Dense[id] > 0 )? tensor_dense18Dense[id] : 0);
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
   std::copy(tensor_dense19bias0bcast, tensor_dense19bias0bcast + 1, tensor_dense19Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense19kernel0, &op_6_ldb, tensor_dense18Relu0, &op_6_lda, &op_6_beta, tensor_dense19Dense, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense19Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense19Dense[id]));
	}
   std::vector<float> ret (tensor_dense19Sigmoid0, tensor_dense19Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
