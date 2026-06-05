//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Fri Jun  5 16:31:09 2026] 

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
std::vector<float> fTensor_dense31bias0 = std::vector<float>(1);
float * tensor_dense31bias0 = fTensor_dense31bias0.data();
std::vector<float> fTensor_dense31kernel0 = std::vector<float>(32);
float * tensor_dense31kernel0 = fTensor_dense31kernel0.data();
std::vector<float> fTensor_dense30bias0 = std::vector<float>(32);
float * tensor_dense30bias0 = fTensor_dense30bias0.data();
std::vector<float> fTensor_dense28kernel0 = std::vector<float>(448);
float * tensor_dense28kernel0 = fTensor_dense28kernel0.data();
std::vector<float> fTensor_dense29kernel0 = std::vector<float>(4096);
float * tensor_dense29kernel0 = fTensor_dense29kernel0.data();
std::vector<float> fTensor_dense29bias0 = std::vector<float>(64);
float * tensor_dense29bias0 = fTensor_dense29bias0.data();
std::vector<float> fTensor_dense28bias0 = std::vector<float>(64);
float * tensor_dense28bias0 = fTensor_dense28bias0.data();
std::vector<float> fTensor_dense30kernel0 = std::vector<float>(2048);
float * tensor_dense30kernel0 = fTensor_dense30kernel0.data();
std::vector<float> fTensor_dense31Sigmoid0 = std::vector<float>(1);
float * tensor_dense31Sigmoid0 = fTensor_dense31Sigmoid0.data();
std::vector<float> fTensor_dense30Dense = std::vector<float>(32);
float * tensor_dense30Dense = fTensor_dense30Dense.data();
std::vector<float> fTensor_dense31bias0bcast = std::vector<float>(1);
float * tensor_dense31bias0bcast = fTensor_dense31bias0bcast.data();
std::vector<float> fTensor_dense30Relu0 = std::vector<float>(32);
float * tensor_dense30Relu0 = fTensor_dense30Relu0.data();
std::vector<float> fTensor_dense31Dense = std::vector<float>(1);
float * tensor_dense31Dense = fTensor_dense31Dense.data();
std::vector<float> fTensor_dense29Relu0 = std::vector<float>(64);
float * tensor_dense29Relu0 = fTensor_dense29Relu0.data();
std::vector<float> fTensor_dense29bias0bcast = std::vector<float>(64);
float * tensor_dense29bias0bcast = fTensor_dense29bias0bcast.data();
std::vector<float> fTensor_dense29Dense = std::vector<float>(64);
float * tensor_dense29Dense = fTensor_dense29Dense.data();
std::vector<float> fTensor_dense28Dense = std::vector<float>(64);
float * tensor_dense28Dense = fTensor_dense28Dense.data();
std::vector<float> fTensor_dense28bias0bcast = std::vector<float>(64);
float * tensor_dense28bias0bcast = fTensor_dense28bias0bcast.data();
std::vector<float> fTensor_dense30bias0bcast = std::vector<float>(32);
float * tensor_dense30bias0bcast = fTensor_dense30bias0bcast.data();
std::vector<float> fTensor_dense28Relu0 = std::vector<float>(64);
float * tensor_dense28Relu0 = fTensor_dense28Relu0.data();


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
   if (tensor_name != "tensor_dense31bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense31bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense31bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense31kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense31kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense31kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense30bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense30bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense30bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense28kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense28kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 448) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 448 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense28kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense29kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense29kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4096) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4096 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense29kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense29bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense29bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense29bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense28bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense28bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 64) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 64 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense28bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense30kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense30kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 2048) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 2048 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense30kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense28bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense28bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense29bias0,{ 64 }, { 1 , 64 });
      std::copy(data, data + 64, tensor_dense29bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense30bias0,{ 32 }, { 1 , 32 });
      std::copy(data, data + 32, tensor_dense30bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense31bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense31bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input8){

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
   std::copy(tensor_dense28bias0bcast, tensor_dense28bias0bcast + 64, tensor_dense28Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense28kernel0, &op_0_ldb, tensor_input8, &op_0_lda, &op_0_beta, tensor_dense28Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense28Relu0[id] = ((tensor_dense28Dense[id] > 0 )? tensor_dense28Dense[id] : 0);
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
   std::copy(tensor_dense29bias0bcast, tensor_dense29bias0bcast + 64, tensor_dense29Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense29kernel0, &op_2_ldb, tensor_dense28Relu0, &op_2_lda, &op_2_beta, tensor_dense29Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 64 ; id++){
      tensor_dense29Relu0[id] = ((tensor_dense29Dense[id] > 0 )? tensor_dense29Dense[id] : 0);
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
   std::copy(tensor_dense30bias0bcast, tensor_dense30bias0bcast + 32, tensor_dense30Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense30kernel0, &op_4_ldb, tensor_dense29Relu0, &op_4_lda, &op_4_beta, tensor_dense30Dense, &op_4_n);

//------ RELU
   for (int id = 0; id < 32 ; id++){
      tensor_dense30Relu0[id] = ((tensor_dense30Dense[id] > 0 )? tensor_dense30Dense[id] : 0);
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
   std::copy(tensor_dense31bias0bcast, tensor_dense31bias0bcast + 1, tensor_dense31Dense);
   BLAS::sgemm_(&op_6_transB, &op_6_transA, &op_6_n, &op_6_m, &op_6_k, &op_6_alpha, tensor_dense31kernel0, &op_6_ldb, tensor_dense30Relu0, &op_6_lda, &op_6_beta, tensor_dense31Dense, &op_6_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense31Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense31Dense[id]));
	}
   std::vector<float> ret (tensor_dense31Sigmoid0, tensor_dense31Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
