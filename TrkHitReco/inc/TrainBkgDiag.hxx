//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Tue Jun  9 20:23:16 2026] 

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
std::vector<float> fTensor_dense23bias0 = std::vector<float>(1);
float * tensor_dense23bias0 = fTensor_dense23bias0.data();
std::vector<float> fTensor_dense23kernel0 = std::vector<float>(4);
float * tensor_dense23kernel0 = fTensor_dense23kernel0.data();
std::vector<float> fTensor_dense22bias0 = std::vector<float>(4);
float * tensor_dense22bias0 = fTensor_dense22bias0.data();
std::vector<float> fTensor_dense22kernel0 = std::vector<float>(32);
float * tensor_dense22kernel0 = fTensor_dense22kernel0.data();
std::vector<float> fTensor_dense21bias0 = std::vector<float>(8);
float * tensor_dense21bias0 = fTensor_dense21bias0.data();
std::vector<float> fTensor_dense21kernel0 = std::vector<float>(24);
float * tensor_dense21kernel0 = fTensor_dense21kernel0.data();
std::vector<float> fTensor_dense23Sigmoid0 = std::vector<float>(1);
float * tensor_dense23Sigmoid0 = fTensor_dense23Sigmoid0.data();
std::vector<float> fTensor_dense23bias0bcast = std::vector<float>(1);
float * tensor_dense23bias0bcast = fTensor_dense23bias0bcast.data();
std::vector<float> fTensor_dense21Relu0 = std::vector<float>(8);
float * tensor_dense21Relu0 = fTensor_dense21Relu0.data();
std::vector<float> fTensor_dense22Relu0 = std::vector<float>(4);
float * tensor_dense22Relu0 = fTensor_dense22Relu0.data();
std::vector<float> fTensor_dense22Dense = std::vector<float>(4);
float * tensor_dense22Dense = fTensor_dense22Dense.data();
std::vector<float> fTensor_dense23Dense = std::vector<float>(1);
float * tensor_dense23Dense = fTensor_dense23Dense.data();
std::vector<float> fTensor_dense21Dense = std::vector<float>(8);
float * tensor_dense21Dense = fTensor_dense21Dense.data();
std::vector<float> fTensor_dense22bias0bcast = std::vector<float>(4);
float * tensor_dense22bias0bcast = fTensor_dense22bias0bcast.data();
std::vector<float> fTensor_dense21bias0bcast = std::vector<float>(8);
float * tensor_dense21bias0bcast = fTensor_dense21bias0bcast.data();


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
   if (tensor_name != "tensor_dense23bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense23bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense23bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense23kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense23kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense23kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense22bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense22bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense22bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense22kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense22kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense22kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense21bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense21bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 8) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 8 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense21bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense21kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense21kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense21kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense21bias0,{ 8 }, { 1 , 8 });
      std::copy(data, data + 8, tensor_dense21bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense22bias0,{ 4 }, { 1 , 4 });
      std::copy(data, data + 4, tensor_dense22bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense23bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense23bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input8){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 1;
   int op_0_n = 8;
   int op_0_k = 3;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 3;
   int op_0_ldb = 8;
   std::copy(tensor_dense21bias0bcast, tensor_dense21bias0bcast + 8, tensor_dense21Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense21kernel0, &op_0_ldb, tensor_input8, &op_0_lda, &op_0_beta, tensor_dense21Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 8 ; id++){
      tensor_dense21Relu0[id] = ((tensor_dense21Dense[id] > 0 )? tensor_dense21Dense[id] : 0);
   }

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 1;
   int op_2_n = 4;
   int op_2_k = 8;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 8;
   int op_2_ldb = 4;
   std::copy(tensor_dense22bias0bcast, tensor_dense22bias0bcast + 4, tensor_dense22Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense22kernel0, &op_2_ldb, tensor_dense21Relu0, &op_2_lda, &op_2_beta, tensor_dense22Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 4 ; id++){
      tensor_dense22Relu0[id] = ((tensor_dense22Dense[id] > 0 )? tensor_dense22Dense[id] : 0);
   }

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 1;
   int op_4_n = 1;
   int op_4_k = 4;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 4;
   int op_4_ldb = 1;
   std::copy(tensor_dense23bias0bcast, tensor_dense23bias0bcast + 1, tensor_dense23Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense23kernel0, &op_4_ldb, tensor_dense22Relu0, &op_4_lda, &op_4_beta, tensor_dense23Dense, &op_4_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense23Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense23Dense[id]));
	}
   std::vector<float> ret (tensor_dense23Sigmoid0, tensor_dense23Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
