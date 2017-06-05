(eval-when (:compile-toplevel :execute :load-toplevel)
  (ql:quickload :cl-cpp-generator))

;; g++ main.cpp -lm -std=gnu++14

;; self testing of fft and other multivariate functions
;; https://ecommons.cornell.edu/bitstream/handle/1813/6243/94-1453.pdf?sequence=1
(in-package :cl-cpp-generator)

(defmacro e (&body body)
  `(statements (<< "std::cout" ,@(loop for e in body collect
				      (cond ((stringp e) `(string ,e))
					    (t e))) "std::endl")))

(let ((code `(with-compilation-unit
		 (include <iostream>)		 
	       (include <cmath>)
	       (include <algorithm>)
	       (include <array>)
	       (include <complex>)

	       ,(let ((n 1024))
		  `(with-compilation-unit
		    (enum Constants (M_MAG_N ,n))

		     (decl (
			   
			    (m_fft_in :type "std::array<std::complex<float>,M_MAG_N>"  :init (list (list ,@ (loop for i below n collect 0.0))))
			    (m_fft_out :type "std::array<std::complex<float>,M_MAG_N>"  :init (list (list ,@ (loop for i below n collect 0.0))))
			   (m_fft_out_mag :type "std::array<float,M_MAG_N>" :init (list (list ,@ (loop for i below n collect 0.0))))
			   ))))

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
	       ;; https://stackoverflow.com/questions/10121574/safe-and-fast-fft
	       (function (fft ((zs :type "std::array<std::complex<float>, N > &"))
			      "template<std::size_t N> void")
			 (let ((j :type "unsigned int" :ctor 0))
			  (dotimes (i (- N 1))
			    (if (< i j)
				(statements
				 (funcall "std::swap" (aref zs i) (aref zs j))))
			    (let ((m :type int :ctor (/ N 2)))
			      (^= j m)
			      (while (== 0 (& j m))
				(/= m 2)
				(^= j m)))))
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

			 (dotimes (i M_MAG_N)
			   (setf (aref m_fft_in i) 0.0
				 (aref m_fft_out i) 0.0
				 (aref m_fft_out_mag i) 0.0))
			 (setf (aref m_fft_in 12) 1.0)
			 (funcall ft m_fft_in m_fft_out)

			 #+nil (dotimes (i M_MAG_N)
			   (setf (aref m_fft_out_mag i) (funcall "std::abs" (aref m_fft_out i))))
			 
			 (dotimes (i M_MAG_N)
			   (macroexpand (e i (string " ") (aref m_fft_out i) (string " ") (funcall "std::exp"
												       (funcall "std::complex<float>"
														0s0
														(/ (* M_PI -2 i)
														   M_MAG_N))))))))))
  (write-source "/home/martin/stage/cl-gen-fft/source/main" "cpp" code))

