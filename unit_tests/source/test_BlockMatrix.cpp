#include <chrono>
#include <functional>
#include <iostream>
#include <iomanip>

#include "simple_math/Utility.h"
#include "simple_math/EigenDefines.h"
#include "simple_math/BlockMatrix.h"
#include "simple_math/BlockVector.h"

using namespace simple_math;

int main()
{
    BlockMatrix<simple_math::GeneralMatrix> block_matrix(2);
    block_matrix.set_block(0, 0, simple_math::GeneralMatrix::Random(2, 2));
    block_matrix.set_block(1, 1, simple_math::GeneralMatrix::Random(2, 2));
    block_matrix.set_block(0, 1, simple_math::GeneralMatrix::Random(2, 2));
    block_matrix.set_block(1, 0, block_matrix.get_block(0, 1).transpose());

    std::cout << block_matrix.str() << std::endl;
	
	return 0;
}