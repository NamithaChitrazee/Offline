//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Tue Jun  9 21:27:09 2026] 

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
std::vector<float> fTensor_dense29bias0 = std::vector<float>(1);
float * tensor_dense29bias0 = fTensor_dense29bias0.data();
std::vector<float> fTensor_dense27kernel0 = std::vector<float>(24);
float * tensor_dense27kernel0 = fTensor_dense27kernel0.data();
std::vector<float> fTensor_dense28kernel0 = std::vector<float>(32);
float * tensor_dense28kernel0 = fTensor_dense28kernel0.data();
std::vector<float> fTensor_dense29kernel0 = std::vector<float>(4);
float * tensor_dense29kernel0 = fTensor_dense29kernel0.data();
std::vector<float> fTensor_dense28bias0 = std::vector<float>(4);
float * tensor_dense28bias0 = fTensor_dense28bias0.data();
std::vector<float> fTensor_dense27bias0 = std::vector<float>(8);
float * tensor_dense27bias0 = fTensor_dense27bias0.data();
std::vector<float> fTensor_dense29Sigmoid0 = std::vector<float>(1);
float * tensor_dense29Sigmoid0 = fTensor_dense29Sigmoid0.data();
std::vector<float> fTensor_dense29Dense = std::vector<float>(1);
float * tensor_dense29Dense = fTensor_dense29Dense.data();
std::vector<float> fTensor_dense28Relu0 = std::vector<float>(4);
float * tensor_dense28Relu0 = fTensor_dense28Relu0.data();
std::vector<float> fTensor_dense28Dense = std::vector<float>(4);
float * tensor_dense28Dense = fTensor_dense28Dense.data();
std::vector<float> fTensor_dense29bias0bcast = std::vector<float>(1);
float * tensor_dense29bias0bcast = fTensor_dense29bias0bcast.data();
std::vector<float> fTensor_dense28bias0bcast = std::vector<float>(4);
float * tensor_dense28bias0bcast = fTensor_dense28bias0bcast.data();
std::vector<float> fTensor_dense27Dense = std::vector<float>(8);
float * tensor_dense27Dense = fTensor_dense27Dense.data();
std::vector<float> fTensor_dense27Relu0 = std::vector<float>(8);
float * tensor_dense27Relu0 = fTensor_dense27Relu0.data();
std::vector<float> fTensor_dense27bias0bcast = std::vector<float>(8);
float * tensor_dense27bias0bcast = fTensor_dense27bias0bcast.data();


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
   if (tensor_name != "tensor_dense29bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense29bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense29bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense27kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense27kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense27kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense28kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense28kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense28kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense29kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense29kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense29kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense28bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense28bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense28bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense27bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense27bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 8) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 8 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense27bias0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense27bias0,{ 8 }, { 1 , 8 });
      std::copy(data, data + 8, tensor_dense27bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense28bias0,{ 4 }, { 1 , 4 });
      std::copy(data, data + 4, tensor_dense28bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense29bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense29bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input10){

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
   std::copy(tensor_dense27bias0bcast, tensor_dense27bias0bcast + 8, tensor_dense27Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense27kernel0, &op_0_ldb, tensor_input10, &op_0_lda, &op_0_beta, tensor_dense27Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 8 ; id++){
      tensor_dense27Relu0[id] = ((tensor_dense27Dense[id] > 0 )? tensor_dense27Dense[id] : 0);
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
   std::copy(tensor_dense28bias0bcast, tensor_dense28bias0bcast + 4, tensor_dense28Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense28kernel0, &op_2_ldb, tensor_dense27Relu0, &op_2_lda, &op_2_beta, tensor_dense28Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 4 ; id++){
      tensor_dense28Relu0[id] = ((tensor_dense28Dense[id] > 0 )? tensor_dense28Dense[id] : 0);
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
   std::copy(tensor_dense29bias0bcast, tensor_dense29bias0bcast + 1, tensor_dense29Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense29kernel0, &op_4_ldb, tensor_dense28Relu0, &op_4_lda, &op_4_beta, tensor_dense29Dense, &op_4_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense29Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense29Dense[id]));
	}
   std::vector<float> ret (tensor_dense29Sigmoid0, tensor_dense29Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
