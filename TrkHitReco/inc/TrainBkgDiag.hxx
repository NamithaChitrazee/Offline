//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Mon May 18 23:55:47 2026] 

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
std::vector<float> fTensor_dense39kernel0 = std::vector<float>(32);
float * tensor_dense39kernel0 = fTensor_dense39kernel0.data();
std::vector<float> fTensor_dense38kernel0 = std::vector<float>(2048);
float * tensor_dense38kernel0 = fTensor_dense38kernel0.data();
std::vector<float> fTensor_dense37bias0 = std::vector<float>(64);
float * tensor_dense37bias0 = fTensor_dense37bias0.data();
std::vector<float> fTensor_dense39bias0 = std::vector<float>(1);
float * tensor_dense39bias0 = fTensor_dense39bias0.data();
std::vector<float> fTensor_dense38bias0 = std::vector<float>(32);
float * tensor_dense38bias0 = fTensor_dense38bias0.data();
std::vector<float> fTensor_dense37kernel0 = std::vector<float>(4096);
float * tensor_dense37kernel0 = fTensor_dense37kernel0.data();
std::vector<float> fTensor_dense36bias0 = std::vector<float>(64);
float * tensor_dense36bias0 = fTensor_dense36bias0.data();
std::vector<float> fTensor_dense36kernel0 = std::vector<float>(448);
float * tensor_dense36kernel0 = fTensor_dense36kernel0.data();
std::vector<float> fTensor_dense39Sigmoid0 = std::vector<float>(1);
float * tensor_dense39Sigmoid0 = fTensor_dense39Sigmoid0.data();
std::vector<float> fTensor_dense39Dense = std::vector<float>(1);
float * tensor_dense39Dense = fTensor_dense39Dense.data();
std::vector<float> fTensor_dense37Relu0 = std::vector<float>(64);
float * tensor_dense37Relu0 = fTensor_dense37Relu0.data();
std::vector<float> fTensor_dense36Dense = std::vector<float>(64);
float * tensor_dense36Dense = fTensor_dense36Dense.data();
std::vector<float> fTensor_dense38Relu0 = std::vector<float>(32);
float * tensor_dense38Relu0 = fTensor_dense38Relu0.data();
std::vector<float> fTensor_dense39bias0bcast = std::vector<float>(1);
float * tensor_dense39bias0bcast = fTensor_dense39bias0bcast.data();
std::vector<float> fTensor_dense38bias0bcast = std::vector<float>(32);
float * tensor_dense38bias0bcast = fTensor_dense38bias0bcast.data();
std::vector<float> fTensor_dense36bias0bcast = std::vector<float>(64);
float * tensor_dense36bias0bcast = fTensor_dense36bias0bcast.data();
std::vector<float> fTensor_dense37Dense = std::vector<float>(64);
float * tensor_dense37Dense = fTensor_dense37Dense.data();
std::vector<float> fTensor_dense37bias0bcast = std::vector<float>(64);
float * tensor_dense37bias0bcast = fTensor_dense37bias0bcast.data();
std::vector<float> fTensor_dense36Relu0 = std::vector<float>(64);
float * tensor_dense36Relu0 = fTensor_dense36Relu0.data();
std::vector<float> fTensor_dense38Dense = std::vector<float>(32);
float * tensor_dense38Dense = fTensor_dense38Dense.data();


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
   if (tensor_name != "tensor_dense39kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense39kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense39kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense38kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense38kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 2048) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 2048 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense38kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense37bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense37bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense37bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense39bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense39bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense39bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense38bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense38bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense38bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense37kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense37kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense37kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense36bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense36bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense36bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense36kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense36kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 448) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 448 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense36kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense36bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense36bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense37bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense37bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense38bias0,{ 32 }, { 1 , 32 });
      std::copy(data, data + 32, tensor_dense38bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense39bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense39bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input10){

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
   std::copy(tensor_dense36bias0bcast, tensor_dense36bias0bcast + 64, tensor_dense36Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense36kernel0, &op_0_ldb, tensor_input10, &op_0_lda, &op_0_beta, tensor_dense36Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense36Relu0[id] = ((tensor_dense36Dense[id] > 0 )? tensor_dense36Dense[id] : 0);
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
   std::copy(tensor_dense37bias0bcast, tensor_dense37bias0bcast + 64, tensor_dense37Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense37kernel0, &op_2_ldb, tensor_dense36Relu0, &op_2_lda, &op_2_beta, tensor_dense37Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense37Relu0[id] = ((tensor_dense37Dense[id] > 0 )? tensor_dense37Dense[id] : 0);
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
   std::copy(tensor_dense38bias0bcast, tensor_dense38bias0bcast + 32, tensor_dense38Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense38kernel0, &op_4_ldb, tensor_dense37Relu0, &op_4_lda, &op_4_beta, tensor_dense38Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 32 ; id++){
      tensor_dense38Relu0[id] = ((tensor_dense38Dense[id] > 0 )? tensor_dense38Dense[id] : 0);
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
   std::copy(tensor_dense39bias0bcast, tensor_dense39bias0bcast + 1, tensor_dense39Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense39kernel0, &op_6_ldb, tensor_dense38Relu0, &op_6_lda, &op_6_beta, tensor_dense39Dense, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense39Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense39Dense[id]));
	}
   std::vector<float> ret (tensor_dense39Sigmoid0, tensor_dense39Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
