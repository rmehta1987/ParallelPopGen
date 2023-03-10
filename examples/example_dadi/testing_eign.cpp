#include <eigen3/Eigen/Dense>
#include <iostream>
 
int main()
{
  Eigen::MatrixXf m(2,2);
  m << 1,2,
       3,4;
 
  std::cout << "1-norm(m)     = " << m.cwiseAbs().colwise().sum().maxCoeff()
            << " == "             << m.colwise().lpNorm<1>().maxCoeff() << std::endl;
 
  std::cout << "infty-norm(m) = " << m.cwiseAbs().rowwise().sum().maxCoeff()
            << " == "             << m.rowwise().lpNorm<1>().maxCoeff() << std::endl;

   Eigen::Vector2f test;
   //calculate mean
   test = m.colwise().mean(); 

   for (int i=0; i < test.size(); i++)
    {
        std::cout << test(i) << std::endl;
    }
    // calculate variance
    Eigen::MatrixXf centered;
    centered = m.rowwise() - test.transpose();
    for (int i=0; i < test.size(); i++)
    {
        for (int j=0; j < test.size(); j++)
    {
        std::cout << centered(i,j) << std::endl;
    }
    }
    //centered = centered.cwiseProduct(m);
     Eigen::Vector2f final;
    final = centered.cwiseProduct(centered).colwise().sum() * 1/2.0;
    

   for (int i=0; i < final.size(); i++)
    {
        std::cout << final(i) << std::endl;
    }


}