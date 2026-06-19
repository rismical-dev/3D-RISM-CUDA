#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include "rism3d.h"

void RISM3D :: cal_Coulomb (string esp) {
  __global__ void coulomb(double * de, double * dfr,
			  double4 * dru, double * dqu,
			  double dx, double dy, double dz,
			  int nx, int ny, int nz, int natu);
  __global__ void fk(double2 *, const double4 * __restrict__ , 
		     const double4 * __restrict__ , const double * __restrict__, 
		     int);
  __global__ void beta(double * dfr, double2 * dfk, double ubeta);
  __global__ void beta2(double * de, double ubeta);

  cout << "synthesizing solute Coulomb potential ..." << endl;
  
  cudaMalloc(&de, ce -> ngrid * sizeof(double));
  cudaMalloc(&dfr, ce -> ngrid * sizeof(double));
  cudaMalloc(&dfk, ce -> ngrid * sizeof(double2));
  cudaMemset(de, 0.0, ce -> ngrid * sizeof(double));
  cudaMemset(dfr, 0.0, ce -> ngrid * sizeof(double));
  cudaMemset(dfk, 0.0, ce -> ngrid * sizeof(double2));

  coulomb <<< g, b >>> (de, dfr, su -> dr, su -> dq,
			ce -> dr[0], ce -> dr[1], ce -> dr[2], 
			ce -> grid[0], ce -> grid[1], ce -> grid[2], su -> num);

  fk <<< g, b >>> (dfk, dgv, su -> dr, su -> dq, su -> num);

  double ubeta = hartree * bohr / (boltzmann * sv -> temper);
  beta <<< g, b >>> (dfr, dfk, ubeta);

  if (esp.empty()) {
    double ubeta = hartree * bohr / (boltzmann * sv -> temper);
    beta2 <<< g, b >>> (de, ubeta);
  } else {

    ifstream in_file(esp.c_str());
    if (!in_file.is_open()) {
      std::cerr << "Error during read esp file: " << esp << std::endl;
      exit(1);
    }

    std::string char80;
    std::getline(in_file, char80);
    std::getline(in_file, char80);

    int natom_cube, di;
    double x0, y0, z0;
    if (!(in_file >> natom_cube >> x0 >> y0 >> z0 >> di)) {
      std::cerr << "Error reading atom count and origin." << std::endl;
      exit(1);
    }
    
    int nx, ny, nz;
    double dx[3], dy[3], dz[3];
    in_file >> nx >> dx[0] >> dx[1] >> dx[2];
    in_file >> ny >> dy[0] >> dy[1] >> dy[2];
    in_file >> nz >> dz[0] >> dz[1] >> dz[2];

    if (nx != ce -> grid[0] || ny != ce -> grid[1] || nz != ce -> grid[2]) {
        std::cerr << "Error. cube file doesn't match input." << std::endl;
	exit(1);
    }

    for (int i = 0; i < natom_cube; ++i) {
      int iatnum;
      double atchg, atx, aty, atz;
      in_file >> iatnum >> atchg >> atx >> aty >> atz;
    }

    double *e = new double[ce -> ngrid];

    for (int ix = 0; ix < ce -> grid[0]; ++ix) {
      for (int iy = 0; iy < ce -> grid[1]; ++iy) {
        for (int iz = 0; iz < ce -> grid[2]; ++iz) {
          double val;
          in_file >> val;
	  size_t k = static_cast<size_t>(ix) 
                   + static_cast<size_t>(iy) * ce -> grid[0]
                   + static_cast<size_t>(iz) * ce -> grid[0] *  ce -> grid[1];
          e[k] = val;
        }
      }
    }
    in_file.close();

    cudaMemcpyAsync(de, e, ce -> ngrid * sizeof(double), cudaMemcpyDefault);
    double ubeta = hartree / (boltzmann * sv -> temper);
    beta2 <<< g, b >>> (de, ubeta);
    delete[] e;
  }
} 


__global__ void coulomb(double * de, double * dfr,
                        double4 * dru, double * dqu,
                        double bx, double by, double bz,
                        int nx, int ny, int nz, int natu) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double rx = ((int)threadIdx.x - nx / 2) * bx;
  double ry = ((int)blockIdx.x - ny / 2) * by;
  double rz = ((int)blockIdx.y - nz / 2) * bz;
  for (int iu = 0; iu < natu; ++iu) {
    double delx = rx - dru[iu].x;
    double dely = ry - dru[iu].y;
    double delz = rz - dru[iu].z;
    double ra = sqrt(delx * delx + dely * dely + delz * delz) ;
    if (ra >= 1.0e-5) {
      double qr = dqu[iu] / ra ;
      de[ip] += qr ;
      dfr[ip] += qr * (1 - exp(- ra)) ;
    } else {
      dfr[ip] += dqu[iu] ;
    }
  }
}


__global__ void fk(double2 * dfk, const double4 * __restrict__ dgv, 
		   const double4 * __restrict__ dru, 
		   const double * __restrict__ dqu, int natu) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  double rk2 = dgv[ip].x * dgv[ip].x
    + dgv[ip].y * dgv[ip].y + dgv[ip].z * dgv[ip].z;
  double rk4i = 1.0 / (rk2 * (rk2 + 1.0));
  for (int iu = 0; iu < natu; ++iu) {
    double ruk = dgv[ip].x * dru[iu].x 
      + dgv[ip].y * dru[iu].y + dgv[ip].z * dru[iu].z;
    double tmp = 4.0 * M_PI * dqu[iu] * rk4i;
    dfk[ip].x += tmp * cos(ruk);
    dfk[ip].y -= tmp * sin(ruk);
  }
}


__global__ void beta(double * dfr, double2 * dfk, double ubeta) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  dfr[ip] *= ubeta;
  dfk[ip].x *= ubeta;
  dfk[ip].y *= ubeta;
}

__global__ void beta2(double * de, double ubeta) {
  unsigned int ip = threadIdx.x + blockIdx.x * blockDim.x
    + blockIdx.y * blockDim.x * gridDim.x;
  de[ip] *= ubeta;
}

