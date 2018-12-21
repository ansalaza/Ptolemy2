package utilities

/**
  * Author: Alex N. Salazar
  * Created on 17-10-2018
  * Contact a.n.salazar@tudelft.nl for questions
  *
  * Description:
  */
object NumericalUtils {

  /**
    * Function to determine the minimum of two numbers
    * @return Boolean
    */
  def min: (Int, Int) => Int = (x, y) => if (x < y) x else y

  /**
    * Function to determine the maximum of two numbers
    * @return Boolean
    */
  def max: (Int, Int) => Int = (x, y) => if (x > y) x else y

}
