//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Tue Oct  6 22:18:37 2026] 

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
std::vector<float> fTensor_dense51bias0 = std::vector<float>(1);
float * tensor_dense51bias0 = fTensor_dense51bias0.data();
std::vector<float> fTensor_dense51kernel0 = std::vector<float>(14);
float * tensor_dense51kernel0 = fTensor_dense51kernel0.data();
std::vector<float> fTensor_dense50bias0 = std::vector<float>(14);
float * tensor_dense50bias0 = fTensor_dense50bias0.data();
std::vector<float> fTensor_dense50kernel0 = std::vector<float>(196);
float * tensor_dense50kernel0 = fTensor_dense50kernel0.data();
std::vector<float> fTensor_dense49bias0 = std::vector<float>(14);
float * tensor_dense49bias0 = fTensor_dense49bias0.data();
std::vector<float> fTensor_dense49kernel0 = std::vector<float>(196);
float * tensor_dense49kernel0 = fTensor_dense49kernel0.data();
std::vector<float> fTensor_dense48bias0 = std::vector<float>(14);
float * tensor_dense48bias0 = fTensor_dense48bias0.data();
std::vector<float> fTensor_dense48kernel0 = std::vector<float>(98);
float * tensor_dense48kernel0 = fTensor_dense48kernel0.data();
std::vector<float> fTensor_dense51bias0bcast = std::vector<float>(1);
float * tensor_dense51bias0bcast = fTensor_dense51bias0bcast.data();
std::vector<float> fTensor_dense51Dense = std::vector<float>(1);
float * tensor_dense51Dense = fTensor_dense51Dense.data();
std::vector<float> fTensor_dense49bias0bcast = std::vector<float>(14);
float * tensor_dense49bias0bcast = fTensor_dense49bias0bcast.data();
std::vector<float> fTensor_dense50Dense = std::vector<float>(14);
float * tensor_dense50Dense = fTensor_dense50Dense.data();
std::vector<float> fTensor_dense50bias0bcast = std::vector<float>(14);
float * tensor_dense50bias0bcast = fTensor_dense50bias0bcast.data();
std::vector<float> fTensor_dense49Dense = std::vector<float>(14);
float * tensor_dense49Dense = fTensor_dense49Dense.data();
std::vector<float> fTensor_dense49Relu0 = std::vector<float>(14);
float * tensor_dense49Relu0 = fTensor_dense49Relu0.data();
std::vector<float> fTensor_dense51Sigmoid0 = std::vector<float>(1);
float * tensor_dense51Sigmoid0 = fTensor_dense51Sigmoid0.data();
std::vector<float> fTensor_dense48Relu0 = std::vector<float>(14);
float * tensor_dense48Relu0 = fTensor_dense48Relu0.data();
std::vector<float> fTensor_dense50Relu0 = std::vector<float>(14);
float * tensor_dense50Relu0 = fTensor_dense50Relu0.data();
std::vector<float> fTensor_dense48Dense = std::vector<float>(14);
float * tensor_dense48Dense = fTensor_dense48Dense.data();
std::vector<float> fTensor_dense48bias0bcast = std::vector<float>(14);
float * tensor_dense48bias0bcast = fTensor_dense48bias0bcast.data();


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
   if (tensor_name != "tensor_dense51bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense51bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense51bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense51kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense51kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 14) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 14 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense51kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense50bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense50bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 14) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 14 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense50bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense50kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense50kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 196) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 196 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense50kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense49bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense49bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 14) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 14 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense49bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense49kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense49kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 196) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 196 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense49kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense48bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense48bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 14) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 14 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense48bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense48kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense48kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 98) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 98 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense48kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense48bias0,{ 14 }, { 1 , 14 });
      std::copy(data, data + 14, tensor_dense48bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense49bias0,{ 14 }, { 1 , 14 });
      std::copy(data, data + 14, tensor_dense49bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense50bias0,{ 14 }, { 1 , 14 });
      std::copy(data, data + 14, tensor_dense50bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense51bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense51bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input13){

//--------- Gemm
   char op_0_transA = 'n';
   char op_0_transB = 'n';
   int op_0_m = 1;
   int op_0_n = 14;
   int op_0_k = 7;
   float op_0_alpha = 1;
   float op_0_beta = 1;
   int op_0_lda = 7;
   int op_0_ldb = 14;
   std::copy(tensor_dense48bias0bcast, tensor_dense48bias0bcast + 14, tensor_dense48Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense48kernel0, &op_0_ldb, tensor_input13, &op_0_lda, &op_0_beta, tensor_dense48Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 14 ; id++){
      tensor_dense48Relu0[id] = ((tensor_dense48Dense[id] > 0 )? tensor_dense48Dense[id] : 0);
   }

//--------- Gemm
   char op_2_transA = 'n';
   char op_2_transB = 'n';
   int op_2_m = 1;
   int op_2_n = 14;
   int op_2_k = 14;
   float op_2_alpha = 1;
   float op_2_beta = 1;
   int op_2_lda = 14;
   int op_2_ldb = 14;
   std::copy(tensor_dense49bias0bcast, tensor_dense49bias0bcast + 14, tensor_dense49Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense49kernel0, &op_2_ldb, tensor_dense48Relu0, &op_2_lda, &op_2_beta, tensor_dense49Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 14 ; id++){
      tensor_dense49Relu0[id] = ((tensor_dense49Dense[id] > 0 )? tensor_dense49Dense[id] : 0);
   }

//--------- Gemm
   char op_4_transA = 'n';
   char op_4_transB = 'n';
   int op_4_m = 1;
   int op_4_n = 14;
   int op_4_k = 14;
   float op_4_alpha = 1;
   float op_4_beta = 1;
   int op_4_lda = 14;
   int op_4_ldb = 14;
   std::copy(tensor_dense50bias0bcast, tensor_dense50bias0bcast + 14, tensor_dense50Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense50kernel0, &op_4_ldb, tensor_dense49Relu0, &op_4_lda, &op_4_beta, tensor_dense50Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 14 ; id++){
      tensor_dense50Relu0[id] = ((tensor_dense50Dense[id] > 0 )? tensor_dense50Dense[id] : 0);
   }

//--------- Gemm
   char op_6_transA = 'n';
   char op_6_transB = 'n';
   int op_6_m = 1;
   int op_6_n = 1;
   int op_6_k = 14;
   float op_6_alpha = 1;
   float op_6_beta = 1;
   int op_6_lda = 14;
   int op_6_ldb = 1;
   std::copy(tensor_dense51bias0bcast, tensor_dense51bias0bcast + 1, tensor_dense51Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense51kernel0, &op_6_ldb, tensor_dense50Relu0, &op_6_lda, &op_6_beta, tensor_dense51Dense, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense51Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense51Dense[id]));
	}
   std::vector<float> ret (tensor_dense51Sigmoid0, tensor_dense51Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
