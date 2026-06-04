import numpy as np
import pyfftw

# This function performs the 3D FFT of a field in single precision
#def FFT3Dr_f(np.array[np.float32_t ndim=3] a, threads):
def FFT3DR_f(a, threads):
	# align arrays
	dims  = len(a)
	a_in  = pyfftw.empty_aligned((dims,dims,dims),    dtype='float32')
	a_out = pyfftw.empty_aligned((dims,dims,dims/2+1),dtype='complex64')

	# plan FFTW
	fftw_plan = pyfftw.FFTW(a_in, a_out, axes=(0,1,2), flags=('FFTW_ESTIMATE',), direction='FFTW_FORWARD', threads=threads)

	# put input array into delta_r and perform FFTW
	a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out


# This function performs the 3D FFT of a field in double precision
def FFT3Dr_d(np.ndarray[np.float64_t,ndim=3] a, int threads):
	# align arrays
	dims  = len(a)
	a_in  = pyfftw.empty_aligned((dims,dims,dims),    dtype='float64')
	a_out = pyfftw.empty_aligned((dims,dims,dims/2+1),dtype='complex128')

	# plan FFTW
	fftw_plan = pyfftw.FFTW(a_in,a_out,axes=(0,1,2), flags=('FFTW_ESTIMATE',), direction='FFTW_FORWARD',threads=threads)

	# put input array into delta_r and perform FFTW
	a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out


# This function performs the 3D FFT of a field in single precision
def IFFT3Dr_f(np.complex64_t[:,:,::1] a, int threads):
	# align arrays
	dims  = len(a)
	a_in  = pyfftw.empty_aligned((dims,dims,dims/2+1),dtype='complex64')
	a_out = pyfftw.empty_aligned((dims,dims,dims),    dtype='float32')

	# plan FFTW
	fftw_plan = pyfftw.FFTW(a_in, a_out, axes=(0,1,2),flags=('FFTW_ESTIMATE',),direction='FFTW_BACKWARD', threads=threads)

	# put input array into delta_r and perform FFTW
	a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out


# This function performs the 3D FFT of a field in double precision
def IFFT3Dr_d(np.complex128_t[:,:,::1] a, int threads):
	# align arrays
	dims  = len(a)
	a_in  = pyfftw.empty_aligned((dims,dims,dims/2+1),dtype='complex128')
	a_out = pyfftw.empty_aligned((dims,dims,dims),    dtype='float64')

	# plan FFTW
	fftw_plan = pyfftw.FFTW(a_in,a_out,axes=(0,1,2),flags=('FFTW_ESTIMATE',),direction='FFTW_BACKWARD',threads=threads)

	# put input array into delta_r and perform FFTW
	a_in [:] = a;  fftw_plan(a_in,a_out);  return a_out
