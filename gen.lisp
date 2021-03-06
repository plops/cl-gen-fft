(eval-when (:compile-toplevel :execute :load-toplevel)
  (ql:quickload :cl-cpp-generator))

;; g++ main.cpp -lm -std=gnu++14

;;
;;https://kfrlib.com/blog/how-c14-and-c17-help-to-write-faster-and-better-code-real-world-examples/

;; self testing of fft and other multivariate functions
;; https://ecommons.cornell.edu/bitstream/handle/1813/6243/94-1453.pdf?sequence=1
(in-package :cl-cpp-generator)

(defmacro e (&body body)
  `(statements (<< "std::cout" ,@(loop for e in body collect
				      (cond ((stringp e) `(string ,e))
					    (t e))) "std::endl")))

(defmacro benchmark (&body body)
  `(let ((start :ctor (funcall current_time)))
     ,@body
     (macroexpand (e "time: " (* (/ 1.0 1000.0) (- (funcall current_time) start)) " ms"))
     ))
(defun rev (x nn)
  (let ((n (floor (log nn 2)))
      (res 0))
  (dotimes (i n)
    (setf (ldb (byte 1 i) res) (ldb (byte 1 (- n 1 i)) x)))
  res))

(let* ((n 32)
      (code `(with-compilation-unit
		 (include <iostream>)
	        (include <iomanip>)		 
	       (include <cmath>)
	       (include <algorithm>)
	       (include <array>)
	       (include <complex>)
	       (include <sys/time.h>) ;; gettimeofday

	       (with-compilation-unit
		     (enum Constants (M_MAG_N ,n))
		 
		   (decl (
			  
			  (m_fft_in :type "std::array<std::complex<float>,M_MAG_N>"  :init (list (list ,@ (loop for i below n collect 0.0))))
			  (m_fft_out :type "std::array<std::complex<float>,M_MAG_N>"  :init (list (list ,@ (loop for i below n collect 0.0))))
			  (m_fft_out2 :type "std::array<std::complex<float>,M_MAG_N>"  :init (list (list ,@ (loop for i below n collect 0.0))))
			  (m_fft_out_mag :type "std::array<float,M_MAG_N>" :init (list (list ,@ (loop for i below n collect 0.0))))
			  )))
	       (function (current_time () "static inline uint64_t")
                                         
                                         (let ((tv :type "struct timeval"))
                                           (funcall gettimeofday &tv nullptr)
                                           (return (+ (* tv.tv_sec 1000000)
                                                      tv.tv_usec))))

	       (function (ft ((in :type "const std::array<std::complex<float>, N > &" )
			      (out :type "std::array<std::complex<float>, N > &" ))
			     "template<std::size_t N> void"
			     )
			 (dotimes (k N)
			   (setf (aref out k) 0))
			 (dotimes (k N)
			   (dotimes (n N)
			     (+= (aref out k) (* (funcall "std::exp"
							  (funcall "std::complex<float>"
								   0s0
								   (/ (* M_PI -2 k n)
								      N)))
						 (aref in n))))))
	       ;; https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
	       (function  (bit_reverse_copy ((in :type "const std::array<std::complex<float>, N > &")
					     (out :type "std::array<std::complex<float>, N > &"))
					    "template<std::size_t N > void")
			  (setf ,@(loop for i below n appending
				  `((aref out ,(rev i n)) (aref in ,i)))))
	       (function
		(fft ((in :type "const std::array<std::complex<float>, N > &")
				(out :type "std::array<std::complex<float>, N > &"))
			       "template<std::size_t N > void")
		(funcall bit_reverse_copy in out)
		,@(loop for s from 1 upto (floor (log n 2)) appending
		       (let ((m (expt 2 s)))
			 `((let ((w_m :type "const auto" :init (funcall "std::complex<float>"
									,(coerce (cos (/ (* pi -2) m)) 'single-float)
									,(coerce (sin (/ (* pi -2) m)) 'single-float))))
			     (for ((k 0) (< k N) (+= k ,m))
				  (let ((w :type "std::complex<float>" :ctor 1))
				    (dotimes (j ,(/ m 2))
				      (let ((t :ctor (* w (aref out (+ k j ,(/ m 2)))))
					    (u :ctor (aref out (+ k j)))
					    )
					(setf (aref out (+ k j)) (+ u t)
					      (aref out (+ k j ,(/ m 2))) (- u t)
					      w (* w w_m)))))))))
		       ))
	       ;; http://users.ece.utexas.edu/~valvano/Starterfiles/FFT.CPP
	       #+nil (function (fft_ ((zs :type "std::array<std::complex<float>, N > &"))
			       "template<std::size_t N > void")
			 (let ((j :ctor 1))
			   (for ((i 1)
				 (< i N)
				 (+= i 2))
				(if (< i j)
				    (statements
				     (funcall ))))))
	       ;; https://stackoverflow.com/questions/10121574/safe-and-fast-fft
	       #+nil
	       (function (fft ((zs :type "std::array<std::complex<float>, N > &"))
			      "template<std::size_t N> void")
			 (let ((j :type "unsigned int" :ctor 0))
			   (raw "// bit reversal")
			   (dotimes (i (- N 1))
			     (if (< i j)
				 (statements
				  (funcall "std::swap" (aref zs i) (aref zs j))))
			     (let ((m :type int :ctor (/ N 2)))
			       (-= j m)
			       (while (== 0 (& j m))
				 (/= m 2)
				 (-= j m)))))
			 (for ((j 1) (< j N) (*= j 2))
			      (dotimes (m j)
				(let ((sign :type "const auto" :ctor 1)
				      (t :type float :ctor (/ (* M_PI sign m) j))
				      (w :ctor #+nil (funcall "std::complex<float>" (funcall cos t) (funcall sin t))
					 (funcall "std::exp" (funcall "std::complex<float>" 0.0 t))))
				  (for ((i m)
					(< i N)
					(+= i (* 2 j)))
				       (let ((zi :ctor (aref zs i))
					     (tz :ctor (* w (funcall zs.at (+ i j)))))
					 (setf (aref zs i) (+ zi (funcall "std::complex<float>" t))
					       (funcall zs.at (+ i j)) (- zi t))))))))
	       
	       (function (main () int)

			 (statements
			  (dotimes (i M_MAG_N)
			    (setf (aref m_fft_in i) 0.0
				  (aref m_fft_out i) 0.0
				  (aref m_fft_out_mag i) 0.0))
			  (setf (aref m_fft_in 1) 1.0)
			  (macroexpand (benchmark
					(dotimes (i 10)
					  (funcall ft m_fft_in m_fft_out)))))

			 (statements
			  (dotimes (i M_MAG_N)
			    (setf (aref m_fft_in i) 0.0
				  (aref m_fft_out2 i) 0.0
				  (aref m_fft_out_mag i) 0.0))
			  (setf (aref m_fft_in 1) 1.0)
			  (macroexpand (benchmark
			     (dotimes (i 10)
			       (funcall fft m_fft_in m_fft_out2)))))
			 
			 #+nil (dotimes  (i M_MAG_N)
				 (setf (aref m_fft_out_mag i) (funcall "std::abs" (aref m_fft_out i))))
			 
			 (dotimes (i M_MAG_N)
			   (macroexpand (e (funcall "std::setw" 6) i (funcall "std::setw" 30) (aref m_fft_out i) (funcall "std::setw" 30) (aref m_fft_out2 i))))))))
  (write-source "/home/martin/stage/cl-gen-fft/source/main" "cpp" code))



