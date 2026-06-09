//Code generated automatically by TMVA for Inference of Model file [TrainBkgDiag.h5] at [Tue Jun  9 17:10:39 2026] 

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
std::vector<float> fTensor_dense11bias0 = std::vector<float>(1);
float * tensor_dense11bias0 = fTensor_dense11bias0.data();
std::vector<float> fTensor_dense10kernel0 = std::vector<float>(32);
float * tensor_dense10kernel0 = fTensor_dense10kernel0.data();
std::vector<float> fTensor_dense10bias0 = std::vector<float>(4);
float * tensor_dense10bias0 = fTensor_dense10bias0.data();
std::vector<float> fTensor_dense11kernel0 = std::vector<float>(4);
float * tensor_dense11kernel0 = fTensor_dense11kernel0.data();
std::vector<float> fTensor_dense9bias0 = std::vector<float>(8);
float * tensor_dense9bias0 = fTensor_dense9bias0.data();
std::vector<float> fTensor_dense9kernel0 = std::vector<float>(24);
float * tensor_dense9kernel0 = fTensor_dense9kernel0.data();
std::vector<float> fTensor_dense11Sigmoid0 = std::vector<float>(1);
float * tensor_dense11Sigmoid0 = fTensor_dense11Sigmoid0.data();
std::vector<float> fTensor_dense11Dense = std::vector<float>(1);
float * tensor_dense11Dense = fTensor_dense11Dense.data();
std::vector<float> fTensor_dense11bias0bcast = std::vector<float>(1);
float * tensor_dense11bias0bcast = fTensor_dense11bias0bcast.data();
std::vector<float> fTensor_dense9Dense = std::vector<float>(8);
float * tensor_dense9Dense = fTensor_dense9Dense.data();
std::vector<float> fTensor_dense10bias0bcast = std::vector<float>(4);
float * tensor_dense10bias0bcast = fTensor_dense10bias0bcast.data();
std::vector<float> fTensor_dense10Relu0 = std::vector<float>(4);
float * tensor_dense10Relu0 = fTensor_dense10Relu0.data();
std::vector<float> fTensor_dense9Relu0 = std::vector<float>(8);
float * tensor_dense9Relu0 = fTensor_dense9Relu0.data();
std::vector<float> fTensor_dense10Dense = std::vector<float>(4);
float * tensor_dense10Dense = fTensor_dense10Dense.data();
std::vector<float> fTensor_dense9bias0bcast = std::vector<float>(8);
float * tensor_dense9bias0bcast = fTensor_dense9bias0bcast.data();


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
   if (tensor_name != "tensor_dense11bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense11bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 1) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 1 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense11bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense10kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense10kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 32) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 32 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense10kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense10bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense10bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense10bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense11kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense11kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 4) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 4 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense11kernel0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense9bias0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense9bias0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 8) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 8 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense9bias0[i];
   f >> tensor_name >> length;
   if (tensor_name != "tensor_dense9kernel0" ) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor name; expected name is tensor_dense9kernel0 , read " + tensor_name;
      throw std::runtime_error(err_msg);
    }
   if (length != 24) {
      std::string err_msg = "TMVA-SOFIE failed to read the correct tensor size; expected size is 24 , read " + std::to_string(length) ;
      throw std::runtime_error(err_msg);
    }
   for (size_t i = 0; i < length; ++i)
      f >> tensor_dense9kernel0[i];
   f.close();
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense9bias0,{ 8 }, { 1 , 8 });
      std::copy(data, data + 8, tensor_dense9bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense10bias0,{ 4 }, { 1 , 4 });
      std::copy(data, data + 4, tensor_dense10bias0bcast);
      delete [] data;
   }
   {
      float * data = TMVA::Experimental::SOFIE::UTILITY::UnidirectionalBroadcast<float>(tensor_dense11bias0,{ 1 }, { 1 , 1 });
      std::copy(data, data + 1, tensor_dense11bias0bcast);
      delete [] data;
   }
}

std::vector<float> infer(float* tensor_input4){

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
   std::copy(tensor_dense9bias0bcast, tensor_dense9bias0bcast + 8, tensor_dense9Dense);
   BLAS::sgemm_(&op_0_transB, &op_0_transA, &op_0_n, &op_0_m, &op_0_k, &op_0_alpha, tensor_dense9kernel0, &op_0_ldb, tensor_input4, &op_0_lda, &op_0_beta, tensor_dense9Dense, &op_0_n);

//------ RELU
   for (int id = 0; id < 8 ; id++){
      tensor_dense9Relu0[id] = ((tensor_dense9Dense[id] > 0 )? tensor_dense9Dense[id] : 0);
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
   std::copy(tensor_dense10bias0bcast, tensor_dense10bias0bcast + 4, tensor_dense10Dense);
   BLAS::sgemm_(&op_2_transB, &op_2_transA, &op_2_n, &op_2_m, &op_2_k, &op_2_alpha, tensor_dense10kernel0, &op_2_ldb, tensor_dense9Relu0, &op_2_lda, &op_2_beta, tensor_dense10Dense, &op_2_n);

//------ RELU
   for (int id = 0; id < 4 ; id++){
      tensor_dense10Relu0[id] = ((tensor_dense10Dense[id] > 0 )? tensor_dense10Dense[id] : 0);
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
   std::copy(tensor_dense11bias0bcast, tensor_dense11bias0bcast + 1, tensor_dense11Dense);
   BLAS::sgemm_(&op_4_transB, &op_4_transA, &op_4_n, &op_4_m, &op_4_k, &op_4_alpha, tensor_dense11kernel0, &op_4_ldb, tensor_dense10Relu0, &op_4_lda, &op_4_beta, tensor_dense11Dense, &op_4_n);
	for (int id = 0; id < 1 ; id++){
		tensor_dense11Sigmoid0[id] = 1 / (1 + std::exp( - tensor_dense11Dense[id]));
	}
   std::vector<float> ret (tensor_dense11Sigmoid0, tensor_dense11Sigmoid0 + 1);
   return ret;
}
};
} //TMVA_SOFIE_TrainBkgDiag

#endif  // ROOT_TMVA_SOFIE_TRAINBKGDIAG
