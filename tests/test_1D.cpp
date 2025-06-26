#include <gtest/gtest.h>
#include "SimpleMesh.hpp"


TEST(SimpleMeshTest, CoordinatesAreCorrect){
    std::vector<double> expected = {0.0, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.0}; 
    SimpleMesh result(0,1,10);
    
    ASSERT_EQ(expected.size(), result.nx());
    for(int i = 0; i < result.nx(); i++){
       EXPECT_DOUBLE_EQ(expected[i], result.x(i)); 
    }
}

TEST(SolutionExact, CorrectExactSolution){
    SimpleMesh mesh(0,1,10);

    std::vector<double> expected = {0.0164611,  0.05475908, 0.10153632, 0.15867018, 0.22845364, 0.31368734, 0.41779202, 0.54494577, 0.7002517,  0.8899428,  1.0}
    std::vector<double> result(12);
    result[11] = 1.0;
    solve_exact(result, &mesh, 0.5, 1.0);

    ASSERT_EQ(expected.size(), result.size())
    for(int i = 0; i < mesh.nx(); i++){
        EXPECT_DOUBLE_EQ(expected[i], result[i]);
    }
}