//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Fri Jun  5 19:28:05 2026] 

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
std::vector<float> fTensor_dense27kernel0 = std::vector<float>(32);
float * tensor_dense27kernel0 = fTensor_dense27kernel0.data();
std::vector<float> fTensor_dense26bias0 = std::vector<float>(32);
float * tensor_dense26bias0 = fTensor_dense26bias0.data();
std::vector<float> fTensor_dense25bias0 = std::vector<float>(64);
float * tensor_dense25bias0 = fTensor_dense25bias0.data();
std::vector<float> fTensor_dense25kernel0 = std::vector<float>(4096);
float * tensor_dense25kernel0 = fTensor_dense25kernel0.data();
std::vector<float> fTensor_dense27bias0 = std::vector<float>(1);
float * tensor_dense27bias0 = fTensor_dense27bias0.data();
std::vector<float> fTensor_dense24bias0 = std::vector<float>(64);
float * tensor_dense24bias0 = fTensor_dense24bias0.data();
std::vector<float> fTensor_dense26kernel0 = std::vector<float>(2048);
float * tensor_dense26kernel0 = fTensor_dense26kernel0.data();
std::vector<float> fTensor_dense24kernel0 = std::vector<float>(640);
float * tensor_dense24kernel0 = fTensor_dense24kernel0.data();
std::vector<float> fTensor_dense27Dense = std::vector<float>(1);
float * tensor_dense27Dense = fTensor_dense27Dense.data();
std::vector<float> fTensor_dense27Sigmoid0 = std::vector<float>(1);
float * tensor_dense27Sigmoid0 = fTensor_dense27Sigmoid0.data();
std::vector<float> fTensor_dense27bias0bcast = std::vector<float>(1);
float * tensor_dense27bias0bcast = fTensor_dense27bias0bcast.data();
std::vector<float> fTensor_dense26Relu0 = std::vector<float>(32);
float * tensor_dense26Relu0 = fTensor_dense26Relu0.data();
std::vector<float> fTensor_dense25Relu0 = std::vector<float>(64);
float * tensor_dense25Relu0 = fTensor_dense25Relu0.data();
std::vector<float> fTensor_dense25Dense = std::vector<float>(64);
float * tensor_dense25Dense = fTensor_dense25Dense.data();
std::vector<float> fTensor_dense26bias0bcast = std::vector<float>(32);
float * tensor_dense26bias0bcast = fTensor_dense26bias0bcast.data();
std::vector<float> fTensor_dense25bias0bcast = std::vector<float>(64);
float * tensor_dense25bias0bcast = fTensor_dense25bias0bcast.data();
std::vector<float> fTensor_dense24Relu0 = std::vector<float>(64);
float * tensor_dense24Relu0 = fTensor_dense24Relu0.data();
std::vector<float> fTensor_dense26Dense = std::vector<float>(32);
float * tensor_dense26Dense = fTensor_dense26Dense.data();
std::vector<float> fTensor_dense24Dense = std::vector<float>(64);
float * tensor_dense24Dense = fTensor_dense24Dense.data();
std::vector<float> fTensor_dense24bias0bcast = std::vector<float>(64);
float * tensor_dense24bias0bcast = fTensor_dense24bias0bcast.data();


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
   if (tensor_name != "tensor_dense27kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense27kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense27kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense26bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense26bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense26bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense25bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense25bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense25bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense25kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense25kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense25kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense27bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense27bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense27bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense24bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense24bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense24bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense26kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense26kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 2048) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 2048 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense26kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense24kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense24kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 640) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 640 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense24kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense24bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense24bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense25bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense25bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense26bias0,{ 32 }, { 1 , 32 });
      std::copy(data, data + 32, tensor_dense26bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense27bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense27bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input7){

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
   std::copy(tensor_dense24bias0bcast, tensor_dense24bias0bcast + 64, tensor_dense24Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense24kernel0, &op_0_ldb, tensor_input7, &op_0_lda, &op_0_beta, tensor_dense24Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense24Relu0[id] = ((tensor_dense24Dense[id] > 0 )? tensor_dense24Dense[id] : 0);
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
   std::copy(tensor_dense25bias0bcast, tensor_dense25bias0bcast + 64, tensor_dense25Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense25kernel0, &op_2_ldb, tensor_dense24Relu0, &op_2_lda, &op_2_beta, tensor_dense25Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense25Relu0[id] = ((tensor_dense25Dense[id] > 0 )? tensor_dense25Dense[id] : 0);
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
   std::copy(tensor_dense26bias0bcast, tensor_dense26bias0bcast + 32, tensor_dense26Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense26kernel0, &op_4_ldb, tensor_dense25Relu0, &op_4_lda, &op_4_beta, tensor_dense26Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 32 ; id++){
      tensor_dense26Relu0[id] = ((tensor_dense26Dense[id] > 0 )? tensor_dense26Dense[id] : 0);
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
   std::copy(tensor_dense27bias0bcast, tensor_dense27bias0bcast + 1, tensor_dense27Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense27kernel0, &op_6_ldb, tensor_dense26Relu0, &op_6_lda, &op_6_beta, tensor_dense27Dense, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense27Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense27Dense[id]));
	}
   std::vector<float> ret (tensor_dense27Sigmoid0, tensor_dense27Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
