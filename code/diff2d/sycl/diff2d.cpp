/****************************************************************
 ****
 **** This file belongs with the course
 **** Parallel Programming in MPI and OpenMP
 **** copyright 2019-2024 Victor Eijkhout eijkhout@tacc.utexas.edu
 ****
 **** diff2d.cpp : 2D diffusion in parallel through SYCL
 ****
 ****************************************************************/

#include <chrono>
#include <iostream>
#include <format>

#include <sycl/sycl.hpp>
using namespace sycl;

int main(int argc,char *argv[])
{
    const int N=atoi(argv[3]),M=atoi(argv[4]);
    int index;
    unsigned long long int ClkPerSec;
    double NSecClk;
    unsigned long int start,end;
    double elapsed_count[10],Average = 0.0;
    float FNorm = 0.0;


    queue cpu_selector(cpu_selector_v);
    queue gpu_selector(gpu_selector_v);
    queue q;
    if(strcmp(argv[2],"cpu") == 0)
        q = cpu_selector;
    else
        q = gpu_selector;

    std::cout << "Device : " << q.get_device().get_info<info::device::name>() << "\n";
    std::cout << "Max Compute Units : " << q.get_device().get_info<info::device::max_compute_units>() << std::endl;

    std::vector<real> Mat_A(N*M,10.0);
    std::vector<real> Mat_Stencil(N*M,0.0);

    buffer<real,2> Buf_a(Mat_A.data(),range<2>(N,M));
    buffer<real,2> Buf_b(Mat_Stencil.data(),range<2>(N,M));
    buffer<real,1> Buf_Fn(&FNorm, range<1>(1));

    Calibrate(&ClkPerSec,NSecClk);

    for(int count = 0;count < atoi(argv[1]);count++)
    {
        start = rdtsc(); 

        // Kernel to compute the 5pt stencil and simultaneously the L2Norm
        q.submit([&] (handler &h)
        {
            accessor D_a(Buf_a,h);
            accessor D_b(Buf_b,h);
            auto D_Fn = reduction(Buf_Fn, h, std::plus<real>());

            h.parallel_for(range<2>(N-2,M-2), D_Fn, [=](item<2> index, auto &sum){
              auto row = index.get_id(0) + 1;
              auto col = index.get_id(1) + 1;

              real stencil_value =
		4*D_a[row][col]
		- D_a[row-1][col] - D_a[row+1][col]
		- D_a[row][col-1] - D_a[row][col+1];
              D_b[row-1][col-1] = stencil_value;
              sum += (stencil_value * stencil_value);
            });
	}).wait();

	FNorm = std::sqrt(FNorm);
	
        q.submit([&] (handler &h)
        {
            accessor D_a(Buf_a,h);
            accessor D_b(Buf_b,h);
            accessor D_Fn(Buf_Fn,h);

            h.parallel_for(range<2>(N-2,M-2), [=](auto index){
              auto row = index.get_id(0) + 1;
              auto col = index.get_id(1) + 1;

              D_a[row][col] = (D_b[row-1][col-1]/D_Fn[0]);
            });
        }).wait();

        end = rdtsc();
    
        elapsed_count[count] = (double)(end - start)/ClkPerSec;
        printf("TTC : %.12f\n",elapsed_count[count]);
    }

    for(int count=1;count<(atoi(argv[1]));count++)
        Average += elapsed_count[count];
   
    std::cout << "\nTime to compute 5pt-Stencil + Power Method (Total) = " << Average << "\n"; 
    Average = Average/(atoi(argv[1]) -1);
    std::cout << "\nTime to compute (Avg over " << atoi(argv[1]) << " loops) = " << Average << "\n";

}
