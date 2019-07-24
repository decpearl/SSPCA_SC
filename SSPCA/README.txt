
*** Structured Sparse Principal Component Analysis (SSPCA) ***  

Matlab package, by Rodolphe Jenatton (www.di.ens.fr/~jenatton/).

For any questions or report any bugs, please contact me at rodolphe [DOT] jenatton [AT] inria [DOT] fr

++DESCRIPTION++

This package is a Matlab program that implements an extension of sparse PCA based on structured sparsity-inducing norms.
When using this toolbox, please reference

    @conference{jenatton2010sspca,
      title={{Structured sparse principal component analysis}},
      author={Jenatton, R. and Obozinski, G. and Bach, F.},
      booktitle={AISTATS},
      year={2010}
    }

For more information, please read the following papers:

    @conference{jenatton2010sspca,
      title={{Structured sparse principal component analysis}},
      author={Jenatton, R. and Obozinski, G. and Bach, F.},
      booktitle={AISTATS},
      year={2010}
    }

    @techreport{jenatton2009,
     title={Structured Variable Selection with Sparsity-Inducing Norms},
     author={Jenatton, R. and Audibert, J.-Y. and Bach, F.},
     institution={arXiv:0904.3523},
     year={2009}
    }



++INSTALLATION++

If need be, specify in SETUP_TOOLBOX.m the path where the files are (if different from the working directory).

At the Matlab prompt, type "SETUP_TOOLBOX" and all mex files will be compiled. 
Detailed demo scripts are:

   * demo_sparse_coding.m
   * demo_group_sparse_coding.m
   * demo_overlapping_groups.m

This version has been tested with R2009a/2010a Matlab version, on 64-bit Linux machines.
Note (see updates below) that some tests have been conducted on Windows 32/64-bits machines: the corresponding files can be found in the directory ModifiedFilesForWin (with the modified SETUP_TOOLBOX.m and the mex files. Note that the mex files should work for both 32/64-bits machines).

The current version is the 1.0 and has been released on February, 20th 2010.

++UPDATES++

%%%% 9/19/2010 %%%%
Bilwaj Gaonkar (Bilwaj.Gaonkar@uphs.upenn.edu) has pointed out the fact that one needs to have blas-devel package installed, not just blas to have it setup under matlab.
Thanks for his notice.

%%%% 11/06/2010 %%%%
Test/Compilation on Windows 32-bits machines by Yahong Han (yahong@zju.edu.cn); thanks for his work!

%%%% 09/02/2011 %%%%
Test/Compilation on Windows 64-bits machines by Shuai Kyle Zheng (kylezheng04@gmail.com); thanks for his contribution! 
His specific setting is windows server 2008 with matlab2010a.

++LICENSE++


This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.


